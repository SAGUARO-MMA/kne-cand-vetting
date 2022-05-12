# JCR May 2022, adapted from Charlie Kilpatrick et al.

from astropy.io import ascii
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sassy_src.models import *
import numpy as np

from sassy_src.models.glade_plus_q3c_orm import GladePlusQ3cRecord
from sassy_src.models.glade_plus_q3c_orm_filters import glade_plus_q3c_orm_filters

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

DB_HOST = os.getenv('DB_HOST', 'localhost')
DB_NAME = os.getenv('DB_NAME', 'sassy')
DB_PASS = os.getenv('DB_PASS', None)
DB_PORT = os.getenv('DB_PORT', 5432)
DB_USER = os.getenv('DB_USER', 'sassy')

DB_CONNECT = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

# Chose radius to be equal to a 100 kpc offset at ~100 Mpc
RADIUS_ARCMIN = 3.183

def galaxy_search(RA: float, Dec: float, _radius: float = RADIUS_ARCMIN, _verbose: bool = False) -> Optional[list]:

    """
    Searches for galaxy match to candidate in GLADE, HECATE, GWGC, SDSS, Legacy DR8 and PS1-STRM catalogs
    Returns 5 galaxies with lowest Pcc values, their offsets, redshifts and magnitudes

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


    for i in range(len(RA)):

        GLADE_matches, GLADE_ra, GLADE_dec, GLADE_offset, GLADE_bmag, GLADE_z, GLADE_z_err = query_GLADE(session, RA[i], Dec[i], _radius)

    _end = time.time()

    print(f"Completed galaxy search in {_end-_begin:.3f} sec")
    print(f"Found {GLADE_matches} of {len(RA)} candidates with a GLADE galaxy match.")

    PCCS = pcc(GLADE_offset,GLADE_bmag)

    best_pcc_args = numpy.where(PCCS == PCCS.min())

    # Write something that sorts all galaxies by Pcc and grabs top 5

    return GLADE_z[best_pcc_args][:5] #gal_name, gal_RA, gal_Dec, pcc, offset, zspec, zphot, zphoterr_low, zphoterr_high, mag, magfilter

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
    sigma = (1/(0.33*np.log(10)))*10**(0.33*(m-24)-2.44)
    prob = 1-np.exp(-(np.pi*(r**2)*sigma))

    return prob

def query_GLADE(session, ra, dec, _radius, _verbose: bool = True):

    """ Query the GLADE+ catalog """
    m=0
    # what if no B-mag?
    gal_offset = []; bmag = []; z = []; z_err = []; gal_ra = []; gal_dec = []

    # set up query
    try:
        query = session.query(GladePlusQ3cRecord)
        query = glade_plus_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}'})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    if len(query.all()) > 0:
        print(len(query.all()))
        m+=1
        for _x in GladePlusQ3cRecord.serialize_list(query.all()):
            gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
            cand = SkyCoord(ra*u.deg, dec*u.deg)
            gal_offset.append(cand.separation(gal).arcsec)
            bmag.append(_x['b'])
            gal_ra.append(_x['ra'])
            gal_dec.append(_x['dec'])
            z.append(_x['z_cmb'])
            z_err.append(_x['z_err'])

    return m, gal_ra, gal_dec, np.array(gal_offset), np.array(bmag), z, z_err

if __name__ == '__main__':

    # get command line argument(s)
    _p = argparse.ArgumentParser(description=f'galaxy_matching.py', formatter_class=argparse.RawTextHelpFormatter)
    _p.add_argument(f'--RA', default=0.0, help="""RA, default '%(default)s'""", type=float, nargs='+')
    _p.add_argument(f'--Dec', default=0.0, help="""DEC., default '%(default)s'""", type=float, nargs='+')
    _p.add_argument(f'--radius', default=RADIUS_ARCMIN, help="""Radius (arcmin), default %(default)s""")
    _p.add_argument(f'--verbose', default=False, action='store_true', help=f'if present, produce more verbose output')
    _a = _p.parse_args()

    # execute
    try:
        output = galaxy_search(RA=_a.RA, Dec =_a.Dec, _radius=float(_a.radius), _verbose=bool(_a.verbose))
        # print(qprob, qso, qoffset)
    except Exception as _x:
        print(f"{_x}")
        print(f"Use:{__doc__}")
