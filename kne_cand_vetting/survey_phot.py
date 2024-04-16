### Query public surveys for photometric detections prior to GW event
### Currently working on ZTF only

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from astropy.time import Time
import os
import requests
from datetime import datetime
import sys
import time
import logging
import numpy as np
import glob
from astropy.io import fits
import argparse
from collections import OrderedDict
import json
import csv
import re

from fundamentals.stats import rolling_window_sigma_clip
from operator import itemgetter
import math

from sassy_q3c_models.ztf_q3c_orm import ZtfQ3cRecord
from sassy_q3c_models.ztf_q3c_orm_filters import ztf_q3c_orm_filters

# constants -- are these still needed?
DB_HOST = os.getenv('POSTGRES_HOST', 'localhost')
DB_NAME = os.getenv('POSTGRES_DB', 'sassy')
DB_PASS = os.getenv('POSTGRES_PASSWORD', None)
DB_PORT = os.getenv('POSTGRES_PORT', 5432)
DB_USER = os.getenv('POSTGRES_USER', 'sassy')

DB_CONNECT = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

RADIUS_ARCSEC = 2.0

logger = logging.getLogger(__name__)

def query_ZTFpubphot(RA: float, Dec: float, _radius: float = RADIUS_ARCSEC, _verbose: bool = False, db_connect: str = DB_CONNECT):
    '''
    Returns public ZTF photometry hosted on sassy.
    '''
    _radius /= 3600 # converting to Degrees
    if _radius <= 0.0:
        raise Exception(f"Invalid input, _radius={_radius}")

    try:
        engine = create_engine(db_connect)
        get_session = sessionmaker(bind=engine)
        session = get_session()
    except Exception as _e2:
        if _verbose:
            print(f"{_e2}")
        raise Exception(f"Failed to connect to database")

    try:
        query = session.query(ZtfQ3cRecord)
        query = ztf_q3c_orm_filters(query, {'cone': f'{RA},{Dec},{_radius}'}) # no time limits at the moment
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({RA}, {Dec})")

    ztfphot = ZtfQ3cRecord.serialize_list(query.all())

    print('{0} photometric detection(s) found in ZTF.'.format(len(ztfphot)))

    return ztfphot

def ATLAS_forcedphot(RA: float, Dec: float, days_ago: float = 200., _verbose: bool = False,
                     token: str = None):
    '''
    Returns stacked and sigma clipped ATLAS forced photometry.
    '''
    BASEURL = "https://fallingstar-data.com/forcedphot"
    # BASEURL = "http://127.0.0.1:8000"

    if token is None:
        token = os.environ.get('ATLAS_API_KEY')
    else:
        print('Using provided token')
    if token is None:
        raise Exception('No token provided')
    else:
        print('Using token from environment')

    headers = {'Authorization': f'Token {token}', 'Accept': 'application/json'}

    t_queryend = Time.now().mjd
    t_querystart = t_queryend - days_ago

    task_url = None
    while not task_url:
        with requests.Session() as s:
            resp = s.post(f"{BASEURL}/queue/", headers=headers, data={'ra': RA, 'dec': Dec,
                                                                      'send_email': False,
                                                                      'mjd_min': t_querystart,
                                                                      'mjd_max': t_queryend,
                                                                      "use_reduced": False,})
            if resp.status_code == 201:  # success
                task_url = resp.json()['url']
                print(f'The task URL is {task_url}')
            elif resp.status_code == 429:  # throttled
                message = resp.json()["detail"]
                print(f'{resp.status_code} {message}')
                t_sec = re.findall(r'available in (\d+) seconds', message)
                t_min = re.findall(r'available in (\d+) minutes', message)
                if t_sec:
                    waittime = int(t_sec[0])
                elif t_min:
                    waittime = int(t_min[0]) * 60
                else:
                    waittime = 10
                print(f'Waiting {waittime} seconds')
                time.sleep(waittime)
            else:
                print(f'ERROR {resp.status_code}')
                print(resp.json())
                sys.exit()

    result_url = None
    taskstarted_printed = False
    while not result_url:
        with requests.Session() as s:
            resp = s.get(task_url, headers=headers)

            if resp.status_code == 200:  # HTTP OK
                if resp.json()['finishtimestamp']:
                    result_url = resp.json()['result_url']
                    print(f"Task is complete with results available at {result_url}")
                elif resp.json()['starttimestamp']:
                    if not taskstarted_printed:
                        print(f"Task is running (started at {resp.json()['starttimestamp']})")
                        taskstarted_printed = True
                    time.sleep(2)
                else:
                    # print(f"Waiting for job to start (queued at {timestamp})")
                    time.sleep(4)
            else:
                print(f'ERROR {resp.status_code}')
                print(resp.text)
                sys.exit()

    with requests.Session() as s:
        textdata = s.get(result_url, headers=headers).text

        # if we'll be making a lot of requests, keep the web queue from being
        # cluttered (and reduce server storage usage) by sending a delete operation
        # s.delete(task_url, headers=headers).json()

    ATLASphot = ATLAS_stack(textdata)

    return ATLASphot

def ATLAS_stack(filecontent, log=logger):
    """
    Function adapted from David Young's :func:`plotter.plot_single_result`
    https://github.com/thespacedoctor/plot-results-from-atlas-force-photometry-service/blob/main/plot_atlas_fp.py
    """
    epochs = ATLAS_read_and_sigma_clip_data(filecontent, log=log)

    # c = cyan, o = arange
    magnitudes = {
        'c': {'mjds': [], 'mags': [], 'magErrs': [], 'lim5sig': []},
        'o': {'mjds': [], 'mags': [], 'magErrs': [], 'lim5sig': []},
        'I': {'mjds': [], 'mags': [], 'magErrs': [], 'lim5sig': []},
    }

    # SPLIT BY FILTER
    for epoch in epochs:
        if epoch["F"] in ["c", "o", "I"]:
            magnitudes[epoch["F"]]["mjds"].append(epoch["MJD"])
            magnitudes[epoch["F"]]["mags"].append(epoch["uJy"])
            magnitudes[epoch["F"]]["magErrs"].append(epoch["duJy"])
            magnitudes[epoch["F"]]['lim5sig'].append(epoch["mag5sig"])

    # STACK PHOTOMETRY IF REQUIRED
    stacked_magnitudes = stack_photometry(magnitudes, binningDays=1)

    return stacked_magnitudes

def ATLAS_read_and_sigma_clip_data(filecontent, log, clippingSigma=2.2):
    """
    Function adapted from David Young's :func:`plotter.read_and_sigma_clip_data`
    https://github.com/thespacedoctor/plot-results-from-atlas-force-photometry-service/blob/main/plot_atlas_fp.py

    *clean up rouge data from the files by performing some basic clipping*
    **Key Arguments:**
    - `fpFile` -- path to single force photometry file
    - `clippingSigma` -- the level at which to clip flux data
    **Return:**
    - `epochs` -- sigma clipped and cleaned epoch data
    """

    # CLEAN UP FILE FOR EASIER READING
    fpData = filecontent.replace("###", "").replace(" ", ",").replace(
        ",,", ",").replace(",,", ",").replace(",,", ",").replace(",,", ",").splitlines()

    # PARSE DATA WITH SOME FIXED CLIPPING
    oepochs = []
    cepochs = []
    csvReader = csv.DictReader(
        fpData, dialect='excel', delimiter=',', quotechar='"')

    for row in csvReader:
        for k, v in row.items():
            try:
                row[k] = float(v)
            except:
                pass
        # REMOVE VERY HIGH ERROR DATA POINTS, POOR CHI SQUARED, OR POOR EPOCHS
        if row["duJy"] > 4000 or row["chi/N"] > 100 or row['mag5sig']<17.:
            continue
        if row["F"] == "c":
            cepochs.append(row)
        if row["F"] == "o":
            oepochs.append(row)

    # SORT BY MJD
    cepochs = sorted(cepochs, key=itemgetter('MJD'), reverse=False)
    oepochs = sorted(oepochs, key=itemgetter('MJD'), reverse=False)

    # SIGMA-CLIP THE DATA WITH A ROLLING WINDOW
    cdataFlux = []
    cdataFlux[:] = [row["uJy"] for row in cepochs]
    odataFlux = []
    odataFlux[:] = [row["uJy"] for row in oepochs]

    maskList = []
    for flux in [cdataFlux, odataFlux]:
        fullMask = rolling_window_sigma_clip(
            log=log,
            array=flux,
            clippingSigma=clippingSigma,
            windowSize=11)
        maskList.append(fullMask)

    try:
        cepochs = [e for e, m in zip(
            cepochs, maskList[0]) if m == False]
    except:
        cepochs = []

    try:
        oepochs = [e for e, m in zip(
            oepochs, maskList[1]) if m == False]
    except:
        oepochs = []

    print('Completed the ``read_and_sigma_clip_data`` function')
    # Returns ordered dictionary of all parameters
    return cepochs + oepochs

def stack_photometry(magnitudes, binningDays=1.):
    """
    Function adapted from David Young's :func:`plotter.stack_photometry`
    https://github.com/thespacedoctor/plot-results-from-atlas-force-photometry-service/blob/main/plot_atlas_fp.py

    *stack the photometry for the given temporal range*
    **Key Arguments:**
        - `magnitudes` -- dictionary of photometry divided into filter sets
        - `binningDays` -- the binning to use (in days)
    **Return:**
        - `summedMagnitudes` -- the stacked photometry
    """

    # IF WE WANT TO 'STACK' THE PHOTOMETRY
    summedMagnitudes = {
        'c': {'mjds': [], 'mags': [], 'magErrs': [], 'n': [], 'lim5sig': []},
        'o': {'mjds': [], 'mags': [], 'magErrs': [], 'n': [], 'lim5sig': []},
        'I': {'mjds': [], 'mags': [], 'magErrs': [], 'n': [], 'lim5sig': []},
    }

    # MAGNITUDES/FLUXES ARE DIVIDED IN UNIQUE FILTER SETS - SO ITERATE OVER
    # FILTERS
    allData = []
    for fil, data in list(magnitudes.items()):
        # WE'RE GOING TO CREATE FURTHER SUBSETS FOR EACH UNQIUE MJD (FLOORED TO AN INTEGER)
        # MAG VARIABLE == FLUX (JUST TO CONFUSE YOU)
        distinctMjds = {}
        for mjd, flx, err, lim in zip(data["mjds"], data["mags"], data["magErrs"], data["lim5sig"]):
            # DICT KEY IS THE UNIQUE INTEGER MJD
            key = str(int(math.floor(mjd / float(binningDays))))
            # FIRST DATA POINT OF THE NIGHTS? CREATE NEW DATA SET
            if key not in distinctMjds:
                distinctMjds[key] = {
                    "mjds": [mjd],
                    "mags": [flx],
                    "magErrs": [err],
                    "lim5sig": [lim]
                }
            # OR NOT THE FIRST? APPEND TO ALREADY CREATED LIST
            else:
                distinctMjds[key]["mjds"].append(mjd)
                distinctMjds[key]["mags"].append(flx)
                distinctMjds[key]["magErrs"].append(err)
                distinctMjds[key]["lim5sig"].append(lim)

        # ALL DATA NOW IN MJD SUBSETS. SO FOR EACH SUBSET (I.E. INDIVIDUAL
        # NIGHTS) ...
        for k, v in list(distinctMjds.items()):
            # GIVE ME THE MEAN MJD
            meanMjd = sum(v["mjds"]) / len(v["mjds"])
            summedMagnitudes[fil]["mjds"].append(meanMjd)
            # GIVE ME THE MEAN FLUX
            meanFLux = sum(v["mags"]) / len(v["mags"])
            summedMagnitudes[fil]["mags"].append(meanFLux)
            # GIVE ME THE COMBINED ERROR
            sum_of_squares = sum(x ** 2 for x in v["magErrs"])
            combError = math.sqrt(sum_of_squares) / len(v["magErrs"])
            summedMagnitudes[fil]["magErrs"].append(combError)
            # 5-sigma limits
            comb5SigLimit = 23.9 - 2.5 * math.log10(5. * combError)
            summedMagnitudes[fil]["lim5sig"].append(comb5SigLimit)
            # GIVE ME NUMBER OF DATA POINTS COMBINED
            n = len(v["mjds"])
            summedMagnitudes[fil]["n"].append(n)
            allData.append({
                'mjd': meanMjd,
                'uJy': meanFLux,
                'duJy': combError,
                'F': fil,
                'n': n,
                'mag5sig': comb5SigLimit
            })
    print('completed the ``stack_photometry`` method')

    return allData

def query_TNSphot(objname: str, BOT_ID: str = None, BOT_NAME: str = None, API_KEY: str = None):
    '''
    Returns photometry from TNS.
    '''

    if BOT_ID is None or BOT_NAME is None or API_KEY is None:
        raise Exception('One or more tokens not provided')
    else:
        print('Using tokens from environment')

    get_obj=[("objname",objname), ("objid",""), ("photometry","1"), ("spectra","0")]
    response=TNS_get(get_obj, BOT_ID, BOT_NAME, API_KEY)

    json_file = is_string_json(response.text)
    json_data=format_to_json(response.text)
    result = json.loads(json_data)

    epochs = result['data']['reply']['photometry']
    allPhot = []
    for e in epochs:
        allPhot.append({
            'jd': e['jd'],
            'mag': e['flux'],
            'magerr': e['fluxerr'],
            'F': e['filters']['name'],
            'limflux': e['limflux'],
            'tel': e['telescope']['name']
        })

    return allPhot

def set_bot_tns_marker(BOT_ID: str = None, BOT_NAME: str = None):
    tns_marker = 'tns_marker{"tns_id": "' + str(BOT_ID) + '", "type": "bot", "name": "' + BOT_NAME + '"}'
    return tns_marker

def TNS_get(get_obj, BOT_ID: str = None, BOT_NAME: str = None, API_KEY: str = None):
    TNS="www.wis-tns.org"
    url_tns_api="https://"+TNS+"/api/get"
    get_url = url_tns_api + "/object"
    tns_marker = set_bot_tns_marker(BOT_ID, BOT_NAME)
    headers = {'User-Agent': tns_marker}
    json_file = OrderedDict(get_obj)
    get_data = {'api_key': API_KEY, 'data': json.dumps(json_file)}
    response = requests.post(get_url, headers = headers, data = get_data)
    return response

def format_to_json(source):                                          #
    # change data to json format and return                          #
    parsed=json.loads(source,object_pairs_hook=OrderedDict)          #
    result=json.dumps(parsed,indent=4)                               #
    return result

def is_string_json(string):
    try:
        json_object = json.loads(string)
    except Exception:
        return False
    return json_object

def SAGUARO_forcedphot(RA: float, Dec: float, CSSField: str):
    '''
    Queries SAGUARO images for photometry of candidates.
    Detections come from SExtractor files and upper limits from headers.
    '''
    radius = 5./3600
    SAGUAROPhot = []
    # Gather individual images from each field --> how to do path to connect sassy to beast?
    files = sorted(glob.glob('http://beast.as.arizona.edu:5013/api/home/data/css/G96/*/*/G96_*'+CSSField+'*sext.gz'))
    for f in files:
        if os.path.exists(f.replace('sext.gz','calb.fz')):
            data_file = f.replace('sext.gz','calb.fz')
        else:
            data_file = f.replace('sext.gz','arch.fz')
        if not os.path.exists(data_file):
    	    continue
        with fits.open(data_file) as hdr:
            header = hdr[1].header
            mjd = header['MJD']
        mag, magerr, flag, ra, dec = np.loadtxt(f,usecols=(3,5,7,11,12),unpack=True)
        target_mag = mag[(ra>RA-radius)&(ra<RA+radius)&(dec>Dec-radius)&(dec<Dec+radius)]
        target_magerr = magerr[(ra>RA-radius)&(ra<RA+radius)&(dec>Dec-radius)&(dec<Dec+radius)]
        target_flag = flag[(ra>RA-radius)&(ra<RA+radius)&(dec>Dec-radius)&(dec<Dec+radius)]
        try:
            SAGUAROPhot.append({'mjd':mjd,
                                'mag':target_mag[0],
                                'magerrs':target_magerr[0],
                                'flags':target_flag[0],
                                })
        except:
            continue

    # Get upper limit from headers, files have been moved to CSS folder on beast
    limitfiles = sorted(glob.glob('http://beast.as.arizona.edu:5013/api/home/data/css/red/*/*/G96_*'+CSSField+'*_trans.fits.fz'))

    # Use most recent observation of field for upper limit
    with fits.open(limitfiles[-1]) as hdr:
        mjd = hdr[1].header['MJD']
        upper_mag = hdr[1].header['T-LMAG']

    SAGUAROPhot.append({'mjd':mjd,
                        'limit':upper_mag,
                        })

    return SAGUAROPhot

if __name__ == '__main__':

    params = argparse.ArgumentParser(description=f'survey_phot.py', formatter_class=argparse.RawTextHelpFormatter)
    params.add_argument(f'--RA', default=0.0, help="""R.A., default '%(default)s'""", type=float)
    params.add_argument(f'--Dec', default=0.0, help="""DEC., default '%(default)s'""", type=float)
    params.add_argument(f'--radius', default=RADIUS_ARCSEC, help="""Radius (arcsec), default %(default)s""")
    # GW event time, if none use current time
    params.add_argument(f'--t_Event', default=datetime.now(), help="""Time of GW event, default %(default)s""" )
    params.add_argument(f'--verbose', default=False, action='store_true', help=f'if present, produce more verbose output')
    a = params.parse_args()

    try:
        # For ZTF query over all time
        ztfphot = query_ZTFpubphot(RA=a.RA, Dec =a.Dec, _radius=float(a.radius), _verbose=bool(a.verbose))
        # For ATLAS + other forced phot, only query -200 days
        atlasphot = ATLAS_forcedphot(RA=a.RA, Dec=a.Dec, t_Event=a.t_Event, _verbose=bool(a.verbose))
    except Exception as _x:
        print(f"{_x}")
        print(f"Use:{__doc__}")
