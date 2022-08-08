# JCR May 2022, adapted from Charlie Kilpatrick, Vic Dong et al.

from astropy.io import ascii
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import numpy

from sassy_q3c_models.glade_plus_q3c_orm import GladePlusQ3cRecord
from sassy_q3c_models.glade_plus_q3c_orm_filters import glade_plus_q3c_orm_filters
from sassy_q3c_models.gwgc_q3c_orm import GwgcQ3cRecord
from sassy_q3c_models.gwgc_q3c_orm_filters import gwgc_q3c_orm_filters
from sassy_q3c_models.hecate_q3c_orm import HecateQ3cRecord
from sassy_q3c_models.hecate_q3c_orm_filters import hecate_q3c_orm_filters
from sassy_q3c_models.sdss12phot_q3c_orm import Sdss12PhotQ3cRecord
from sassy_q3c_models.sdss12phot_q3c_orm_filters import sdss12phot_q3c_orm_filters

from typing import Optional
from astropy.coordinates import SkyCoord
from astropy import cosmology
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import z_at_value
cosmo = FlatLambdaCDM(H0=69.6 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

import argparse
import os
import time

# +
# __doc__
# -
__doc__ = """ PYTHONPATH=/home/phil_daly/SASSyII python3 catalogs.py --help """

DB_HOST = os.getenv('POSTGRES_HOST', 'localhost')
DB_NAME = os.getenv('POSTGRES_DB', 'sassy')
DB_PASS = os.getenv('POSTGRES_PASSWORD', None)
DB_PORT = os.getenv('POSTGRES_PORT', 5432)
DB_USER = os.getenv('POSTGRES_USER', 'sassy')

DB_CONNECT = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

# Chose radius to be equal to a 100 kpc offset at ~100 Mpc
RADIUS_ARCMIN = 3.183
PCC_THRESHOLD = 0.80

def galaxy_search(RA: float, Dec: float, _radius: float = RADIUS_ARCMIN, _pcc_thresh: float = PCC_THRESHOLD, _verbose: bool = False) -> Optional[list]:

    """
    Searches for galaxy match to candidate in GLADE, HECATE, GWGC, SDSS, Legacy DR8 and PS1-STRM catalogs
    Returns 10 galaxies with lowest Pcc values (must be < PCC_THRESHOLD), their offsets, redshifts and magnitudes

    PARAMETERS
    ----------
    RA, Dec : transient coordinates, list of floats
        degrees
    _radius : search radius in arcminutes

    RETURNS
    -------
    gal_name : catalogued name of galaxy, list of str
    gal_RA, gal_Dec : RA and Dec of galaxies, list of floats
        degrees
    pcc : probability of chance coincidence, list of floats
    offset : transient-galaxy offset, list of floats
        arcseconds
    zspec : spectroscopic redshift, -99.0 if none, list of floats
    zphot : photometric redshift, -99.0 if none, list of floats
    etc.
    """

    _radius /= 60.0 # convert to degrees

    # connect to database
    try:
        engine = create_engine(DB_CONNECT)
        get_session = sessionmaker(bind=engine)
        session = get_session()
    except Exception as _e2:
        if _verbose:
            print(f"{_e2}")
        raise Exception(f"Failed to connect to database")

    # loop through RA, Dec here
    _begin = time.time()
    glade=0; gwgc=0; hecate=0; sdss=0

    for i in range(len(RA)):

        # Find matches in GLADE:
        GLADE_matches, GLADE_ra, GLADE_dec, GLADE_offset, GLADE_mag, GLADE_filt, GLADE_dist, GLADE_dist_err, GLADE_distflag, GLADE_source, GLADE_name = query_GLADE(session, RA[i], Dec[i], _radius)
        if GLADE_matches>0:
            glade+=1

        # Find matches in GWGC:
        GWGC_matches, GWGC_ra, GWGC_dec, GWGC_offset, GWGC_mag, GWGC_filt, GWGC_dist, GWGC_dist_err, GWGC_source, GWGC_name = query_GWGC(session, RA[i], Dec[i], _radius)
        if GWGC_matches>0:
            gwgc+=1

        HECATE_matches, HECATE_ra, HECATE_dec, HECATE_offset, HECATE_mag, HECATE_filt, HECATE_dist, HECATE_dist_err, HECATE_source, HECATE_name = query_hecate(session, RA[i], Dec[i], _radius)
        if HECATE_matches>0:
            hecate+=1

        SDSS_matches, SDSS_ra, SDSS_dec, SDSS_offset, SDSS_mag, SDSS_filt, SDSS_dist, SDSS_dist_err, SDSS_source = query_sdss12phot(session, RA[i], Dec[i], _radius)
        if SDSS_matches>0:
            sdss+=1

        # sum the findings, turn into numpy arrays
        tot_names = numpy.array(GLADE_name + GWGC_name + HECATE_name)
        tot_offsets = numpy.array(GLADE_offset + GWGC_offset + HECATE_offset)
        tot_mags = numpy.array(GLADE_mag + GWGC_mag + HECATE_mag)
        tot_ra = numpy.array(GLADE_ra + GWGC_ra + HECATE_ra)
        tot_dec = numpy.array(GLADE_dec + GWGC_dec + HECATE_dec)
        tot_filt = numpy.array(GLADE_filt + GWGC_filt + HECATE_filt)
        tot_dists = numpy.array(GLADE_dist + GWGC_dist + HECATE_dist)
        tot_dist_errs = numpy.array(GLADE_dist_err + GWGC_dist_err + HECATE_dist_err)
        tot_source = numpy.array(GLADE_source + GWGC_source + HECATE_source)

        PCCS = pcc(tot_offsets,tot_mags)

        # put some basic cut on Pcc ?
        pcc_args = numpy.argsort(PCCS)[:10]
        cond = (PCCS[pcc_args] < _pcc_thresh)

    _end = time.time()

    print(f"Completed galaxy search in {_end-_begin:.3f} sec")
    print(f"Found {glade} of {len(RA)} candidates with a GLADE galaxy match.")
    print(f"Found {gwgc} of {len(RA)} candidates with a GWGC galaxy match.")
    print(f"Found {hecate} of {len(RA)} candidates with a HECATE galaxy match.")
    print(f"Found {sdss} of {len(RA)} candidates with a SDSS DR12 Photo-z Catalog galaxy match.")

    all_data = [{'ID':tot_names[pcc_args][cond][i],'PCC':PCCS[pcc_args][cond][i],'RA':tot_ra[pcc_args][cond][i],'Dec':tot_dec[pcc_args][cond][i],'Dist':tot_dists[pcc_args][cond][i],'DistErr':tot_dist_errs[pcc_args][cond][i],'Mags':tot_mags[pcc_args][cond][i],'Filter':tot_filt[pcc_args][cond][i],'Source':tot_source[pcc_args][cond][i]} for i in range(len(PCCS[pcc_args][cond]))]

    return len(PCCS[pcc_args][cond]), all_data

def pcc(r,m):
    """
    Probability of chance coincidence calculation (Bloom et al. 2002)

    PARAMETERS
    ----------
    r : transient-galaxy offsets, array of floats
        arcseconds
    m : magnitudes of galaxies, array of floats

    RETURNS
    -------
    Pcc values : array of floats [0,1]
    """
    sigma = (1/(0.33*numpy.log(10)))*10**(0.33*(m-24)-2.44)
    prob = 1-numpy.exp(-(numpy.pi*(r**2)*sigma))

    return prob

def sort_names(catalog,_dict):
    """Chooses preferred name to display"""

    if catalog=='GLADE':
        keys = ['pgc','hyperleda','gwgc','wise','twomass','sdss','gid']
    elif catalog=='GWGC':
        keys = ['name', 'pgc', 'gid']
    elif catalog=='HECATE':
        keys = ['objname', 'id_ned', 'id_iras', 'pgc', 'id_2mass', 'sdss_specid', 'sdss_photid', 'hid']
    else:
        print('No catalog name keys...')

    for key in keys:
        if _dict[key]!='null' and _dict[key]!=-1 and _dict[key]!='-':
            if key=='pgc':
                name = 'PGC'+str(_dict[key])
            elif key=='wise' or key=='twomass':
                name = key+_dict[key]
            else:
                name = _dict[key]
            break
        if key==keys[-1]:
            name = catalog+str(_dict[key])

    return name

def query_GLADE(session, ra, dec, _radius, _verbose: bool = True):

    """
    Query the GLADE+ catalog
    ----------------
    DIST FLAG PARAMETERS:
    0 = the galaxy has no measured redshift or distance value
    1 = it has a measured photometric redshift from which we have calculated its luminosity distance
    2: it has a measured luminosity distance value from which we have calculated its redshift
    3: it has a measured spectroscopic redshift from which we have calculated its luminosity distance
    """
    m=0
    # what if no B-mag?
    gal_offset = []; mag = []; filt = []; dist = []; dist_err = []; gal_ra = []; gal_dec = []; distflag = []; source = []; name = []

    # set up query
    try:
        query = session.query(GladePlusQ3cRecord)
        query = glade_plus_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}'})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    if len(query.all()) > 0:
        m+=1
        for _x in GladePlusQ3cRecord.serialize_list(query.all()):
            if _x['b']== _x['b']:
                mag.append(_x['b'])
                filt.append('B')
                gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['ra']) # degrees
                gal_dec.append(_x['dec']) # degrees
                dist.append(_x['d_l']) # Mpc
                dist_err.append(_x['d_l_err']) # Mpc
                distflag.append(_x['dist_flag'])
                source.append('GLADE')
                name.append(sort_names('GLADE',_x))
            elif _x['b_j']== _x['b_j']:
                mag.append(_x['b_j'])
                filt.append('B_j')
                gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['ra'])
                gal_dec.append(_x['dec'])
                dist.append(_x['d_l'])
                dist_err.append(_x['d_l_err'])
                distflag.append(_x['dist_flag'])
                source.append('GLADE')
                best_id, best_id_source = sort_names('GLADE',_x)
                name.append(best_id)
                name.append(sort_names('GLADE',_x))
            # else:
            #     mag.append(_x['w1'])
            #     filt.append('W1')
            #     gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
            #     cand = SkyCoord(ra*u.deg, dec*u.deg)
            #     gal_offset.append(cand.separation(gal).arcsec)
            #     gal_ra.append(_x['ra'])
            #     gal_dec.append(_x['dec'])
            #     z.append(_x['z_cmb'])
            #     z_err.append(_x['z_err'])

    return m, gal_ra, gal_dec, gal_offset, mag, filt, dist, dist_err, distflag, source, name

def query_GWGC(session, ra, dec, _radius, _verbose: bool = True):

    """ Query the GWGC catalog """
    m=0
    # what if no B-mag?
    gal_offset = []; mag = []; filt = []; dist = []; dist_err = []; gal_ra = []; gal_dec = []; distflag = []; source = []
    name = []

    # set up query
    try:
        query = session.query(GwgcQ3cRecord)
        query = gwgc_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}'})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    if len(query.all()) > 0:
        m+=1
        for _x in GwgcQ3cRecord.serialize_list(query.all()):
            if _x['b_app']== _x['b_app']:
                mag.append(_x['b_app'])
                filt.append('B')
                gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['ra'])
                gal_dec.append(_x['dec'])
                dist.append(_x['dist']) # Mpc
                dist_err.append(_x['e_dist']) # Mpc
                source.append('GWGC')
                name.append(sort_names('GWGC',_x))

    return m, gal_ra, gal_dec, gal_offset, mag, filt, dist, dist_err, source, name

def query_hecate(session, ra, dec, _radius, _verbose: bool = True):

    # Query the HECATE catalog
    m=0
    gal_offset = []; mag = []; filt = []; dist = []; dist_err = []; gal_ra = []; gal_dec = []; distflag = []; source = []
    name = []

    # set up query
    try:
        query = session.query(HecateQ3cRecord)
        query = hecate_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}'})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    if len(query.all()) > 0:
        m+=1
        for _x in HecateQ3cRecord.serialize_list(query.all()):
            if _x['bt']== _x['bt']:
                mag.append(_x['bt'])
                filt.append('B')
                gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['ra'])
                gal_dec.append(_x['dec'])
                dist.append(_x['d']) # Mpc
                dist_err.append(_x['e_d']) # Mpc
                source.append('HECATE')
                name.append(sort_names('HECATE',_x))

    return m, gal_ra, gal_dec, gal_offset, mag, filt, dist, dist_err, source, name

def query_sdss12phot(session, ra, dec, _radius, _verbose: bool = True):

    # Query the SDSS DR12 Photo-z Catalog

    m=0
    gal_offset = []; mag = []; filt = []; dist = []; dist_err = []; gal_ra = [];
    gal_dec = []; distflag = []; source = []; dist_q = []

    try:
        query = session.query(Sdss12PhotQ3cRecord)
        query = sdss12phot_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}'})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    # print(query.all())
    if len(query.all()) > 0:
        m+=1
        print('SDSS, ', Sdss12PhotQ3cRecord.serialize_list(query.all())[0])
        for _x in Sdss12PhotQ3cRecord.serialize_list(query.all()):
            if _x['rmag']== _x['rmag']:
                mag.append(_x['rmag'])
                filt.append('r')
                gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['ra'])
                gal_dec.append(_x['dec'])

                # print(_x)

                # Split spec and photo-z:
                # if _x['z_sp']
                # dist.append(_x['d']) # Mpc
                # dist_err.append(_x['e_d']) # Mpc
                # source.append('SDSS_DR12')

    return m, gal_ra, gal_dec, gal_offset, mag, filt, dist, dist_err, source


if __name__ == '__main__':

    # get command line argument(s)
    _p = argparse.ArgumentParser(description=f'galaxy_matching.py', formatter_class=argparse.RawTextHelpFormatter)
    _p.add_argument(f'--RA', default=0.0, help="""RA, default '%(default)s'""", type=float, nargs='+')
    _p.add_argument(f'--Dec', default=0.0, help="""DEC., default '%(default)s'""", type=float, nargs='+')
    _p.add_argument(f'--radius', default=RADIUS_ARCMIN, help="""Radius (arcmin), default %(default)s""")
    _p.add_argument(f'--pcc_threshold', default=PCC_THRESHOLD, help="""Maximum Probability of Chance Coincidence Value to Return Galaxies, default %(default)s""")
    _p.add_argument(f'--verbose', default=False, action='store_true', help=f'if present, produce more verbose output')
    _a = _p.parse_args()

    # execute
    try:
        output = galaxy_search(RA=_a.RA, Dec =_a.Dec, _radius=float(_a.radius), _pcc_thresh=float(_a.pcc_threshold), _verbose=bool(_a.verbose))
        # print(qprob, qso, qoffset)
    except Exception as _x:
        print(f"{_x}")
        print(f"Use:{__doc__}")
