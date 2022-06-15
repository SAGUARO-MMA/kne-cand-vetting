#!/usr/bin/env python3


# +
# import(s)
# -
from astropy.io import ascii
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sassy_q3c_models import *
from sassy_q3c_models.milliquas_q3c_orm import MilliQuasQ3cRecord
from sassy_q3c_models.milliquas_q3c_orm_filters import milliquas_q3c_orm_filters
from sassy_q3c_models.asassn_q3c_orm import AsAssnQ3cRecord
from sassy_q3c_models.asassn_q3c_orm_filters import asassn_q3c_orm_filters

from typing import Optional
from astropy.coordinates import SkyCoord
from astropy import units as u

import argparse
import os
import time


# +
# __doc__
# -
__doc__ = """ PYTHONPATH=/home/phil_daly/SASSyII python3 catalogs.py --help """


# +
# constant(s)
# -
DB_HOST = os.getenv('DB_HOST', 'localhost')
DB_NAME = os.getenv('DB_NAME', 'sassy')
DB_PASS = os.getenv('POSTGRES_PASSWORD', None)
DB_PORT = os.getenv('DB_PORT', 5432)
DB_USER = os.getenv('DB_USER', 'sassy')

DB_CONNECT = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
FILE = os.path.abspath(os.path.expanduser("comb_master_FINAL.dat"))
RADIUS_ARCSEC = 2.0

# +
# function: static_cats_query()
# -
#  ->  allows you to optionall accept inputs as lists
def static_cats_query(RA: float, Dec: float, _radius: float = RADIUS_ARCSEC, _verbose: bool = False) -> Optional[list]:

    # check input(s)
    # _file = os.path.abspath(os.path.expanduser(f"{_file}"))
    # if not os.path.exists(_file):
    #     raise Exception(f"invalid input, _file={_file}")

    _radius /= 3600.0 # convert to degrees
    if _radius <= 0.0:
        raise Exception(f"invalid input, _radius={_radius}")

    # read the file into an astropy ascii table
    # _candidates = None
    # try:
    #     _candidates = ascii.read(_file)
    # except Exception as _e1:
    #     if _verbose:
    #         print(f"{_e1}")
    #     raise Exception(f"Unable to read {_file}")
    # else:
    #     if len(_candidates) <= 0:
    #         raise Exception(f"File {_file} contains no data")

    # connect to database
    try:
        engine = create_engine(DB_CONNECT)
        get_session = sessionmaker(bind=engine)
        session = get_session()
    except Exception as _e2:
        if _verbose:
            print(f"{_e2}")
        raise Exception(f"Failed to connect to database")

    # _names = [_ for _ in _candidates['Name']]

    # gather co-ordinates into (Ra, Dec) tuples
    _coords = tuple(zip(RA, Dec))
    _names = 'default'

    # Pass to milliquas_query
    qprob, qso, qoffset = milliquas_query(session, _coords, _names, _radius)

    asassnprob, asassn, asassnoffset = asassn_query(session, _coords, _names, _radius)

    session.close()

    return qprob, qso, qoffset, asassnprob, asassn, asassnoffset

def milliquas_query(session, coords, names, _radius, _verbose: bool = True):
    """ Query the Million Quasar Catalog (Flesch 2021) for matches to kilonova candidates """

    # set variable(s)
    _begin = time.time()

    qso = []; qoffset = []; qprob = []
    match=0

    # for all (RA, Dec) tuples, execute cone search and log candidates with qpct > 97(%)
    for _i, _e in enumerate(coords):

        # set up query
        try:
            query = session.query(MilliQuasQ3cRecord)
            query = milliquas_q3c_orm_filters(query, {'cone': f'{_e[0]},{_e[1]},{_radius}', 'q__gte': 97})
        except Exception as _e3:
            if _verbose:
                print(f"{_e3}")
            print(f"Failed to execute query for RA, Dec = ({_e[0]}, {_e[1]}), index={_i}")
            continue

        # execute the query
        if len(query.all()) == 0:
            qso.append(None)
            qprob.append(0.0)
            qoffset.append(None)
        else:
            match+=1
            for _x in MilliQuasQ3cRecord.serialize_list(query.all()):
                print(f'>>> QUASAR MATCH at RA, Dec = ({_e[0]},{_e[1]}), index={_i}!')

                # add the query dictionary to the qso list and modify the qprob list
                qso.append(_x['name']) #{**_x, **{'Candidate': names[_i], 'Probability': 1.0, 'Candidate_RA': _e[0], 'Candidate_Dec': _e[1]}})
                qprob.append(1.0)

                QSO = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(_e[0]*u.deg, _e[1]*u.deg)
                qoffset.append(cand.separation(QSO).arcsec)

    # done
    _end = time.time()

    # return list of probabilities (although I think you'd be better off returning the qso list!)
    # print(f"QSO Candidates: {qso}")
    # print(f"QSO Probabilities: {qprob}")
    print(f"Completed Milliquas search in {_end-_begin:.3f} sec")
    print(f"Found {match} QSOs in {len(coords)} candidates")

    return qprob, qso, qoffset

def asassn_query(session, coords, names, _radius, _verbose: bool = False):
    """ Query the ASAS-SN variable star catalog (Flesch 2021) for matches to kilonova candidates """

    _begin = time.time()

    starprob = []; star = []; staroffset = []
    match=0

    for _i, _e in enumerate(coords):

        # set up query
        try:
            query = session.query(AsAssnQ3cRecord)
            query = asassn_q3c_orm_filters(query, {'cone': f'{_e[0]},{_e[1]},{_radius}'})
        except Exception as _e3:
            if _verbose:
                print(f"{_e3}")
            print(f"Failed to execute query for RA, Dec = ({_e[0]}, {_e[1]}), index={_i}")
            continue
        # execute the query
        if len(query.all()) == 0:
            star.append(None)
            starprob.append(0.0)
            staroffset.append(None)
        else:
            match+=1
            for _x in AsAssnQ3cRecord.serialize_list(query.all()):
                print(f'>>> ASAS-SN Variable Star MATCH at RA, Dec = ({_e[0]},{_e[1]}), index={_i}!')

                # add the query dictionary to the qso list and modify the qprob list
                star.append(_x['asassn_name']) #{**_x, **{'Candidate': names[_i], 'Probability': 1.0, 'Candidate_RA': _e[0], 'Candidate_Dec': _e[1]}})
                starprob.append(1.0)

                asassn = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(_e[0]*u.deg, _e[1]*u.deg)
                staroffset.append(cand.separation(asassn).arcsec)

    _end = time.time()

    print(f"Completed ASAS-SN search in {_end-_begin:.3f} sec")
    print(f"Found {match} variable stars in {len(coords)} candidates")

    return starprob, star, staroffset


# +
# main()
# -
if __name__ == '__main__':

    # get command line argument(s)
    _p = argparse.ArgumentParser(description=f'catalogs.py', formatter_class=argparse.RawTextHelpFormatter)
    _p.add_argument(f'--RA', default=0.0, help="""RA, default '%(default)s'""", type=float, nargs='+')
    _p.add_argument(f'--Dec', default=0.0, help="""DEC., default '%(default)s'""", type=float, nargs='+')
    _p.add_argument(f'--radius', default=RADIUS_ARCSEC, help="""Radius (arcsec), default %(default)s""")
    _p.add_argument(f'--verbose', default=False, action='store_true', help=f'if present, produce more verbose output')
    _a = _p.parse_args()

    # execute
    try:
        asassnprob, asassn, asassnoffset, qprob, qso, qoffset = static_cats_query(RA=_a.RA, Dec =_a.Dec, _radius=float(_a.radius), _verbose=bool(_a.verbose))
        # print(qprob, qso, qoffset)
    except Exception as _x:
        print(f"{_x}")
        print(f"Use:{__doc__}")



### DOES NOT WORK
def gaia_query(ra,dec,rad=2):
    '''
    Query Gaia eDR3 for a variable star match

    PARAMETERS
    ----------
    ra, dec : array of floats
        degrees
    radius : float
        default = 2, arcseconds

    RETURNS
    -------
    match : array of floats
        1. = match to Gaia source meeting criteria 0. = no match

    '''
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    print(coords)

    # Query statement --> recall you're using an array

    dist = (r['dist'].data*u.deg).to(u.arcsec) # perhaps distance from center? units = degrees?
    pmra = r['pmra'] # proper motion, RA in mas/yr
    pmdec = r['pmdec'] # proper motion, Dec in mas/yr
    pmra_err = r['pmra_error'] # in mas/yr
    pmdec_err = r['pmdec_error'] # in mas/yr
    variable_flag = r['phot_variable_flag'] # flag for variable stars
    parallax_sig = r['parallax_over_error'] # parallax over error

    if len(dist) != 0: # meaning, a match is found

        pmra=pmra[0]
        pmdec=pmdec[0]
        pmra_err=pmra_err[0]
        pmdec_err=pmdec_err[0]
        parallax_sig = parallax_sig[0]

        pm_tot = (pmra**2 + pmdec**2)**.5 # find magnitude of proper motion

        err_pm_tot = (pmra_err**2 + pmdec_err**2)**.5 # find magnitude of proper motion error

        vf = variable_flag[0].decode("utf-8") # make variable flag readable

        # conditions for ruling a match out as a star:
        # 1) proper motion > 3 x PM error
        # 2) parallax significance > 8 (more details in Tachibana & Miller 2018)
        # 3) Gaia has flagged source as a variable star (very few)

        if pm_tot > 3*err_pm_tot or parallax_sig > 8 or vf=='VARIABLE':
            print('star: ',vf,pm_tot,err_pm_tot)
            return 1.0

        else:
            return 0.0 # match in Gaia, but did not meet criteria of being stellar

    else:
        return -99.0 # no match in Gaia

    return coords # q
