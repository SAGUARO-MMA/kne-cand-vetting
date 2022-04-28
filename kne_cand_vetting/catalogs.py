#!/usr/bin/env python3


# +
# import(s)
# -
from astropy.io import ascii
from sassy_src.models.milliquas_q3c_orm import *
from sassy_src.models.milliquas_q3c_orm_cli import *
from sassy_src.models.milliquas_q3c_orm_filters import *
# from sassy_src.models.asassn_q3c_orm import *
# from sassy_src.models.asassn_q3c_orm_filters import *
# from sassy_src.models.asassn_q3c_orm_cli import *
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
DB_CONNECT = f"postgresql+psycopg2://sassy:SASSy_520@localhost:5432/sassy"
FILE = os.path.abspath(os.path.expanduser("comb_master_FINAL.dat"))
RADIUS_ARCSEC = 2.0


# +
# function: milliquas_query()
# -
#  ->  allows you to optionall accept inputs as lists
def static_cats_query(_file: str = FILE, _radius: float = RADIUS_ARCSEC, _verbose: bool = False) -> Optional[list]:
    """ Query the Million Quasar Catalog (Flesch 2021) for matches to kilonova candidates """

    # check input(s)
    _file = os.path.abspath(os.path.expanduser(f"{_file}"))
    if not os.path.exists(_file):
        raise Exception(f"invalid input, _file={_file}")

    _radius /= 3600.0 # convert to degrees
    if _radius <= 0.0:
        raise Exception(f"invalid input, _radius={_radius}")

    # read the file into an astropy ascii table
    _candidates = None
    try:
        _candidates = ascii.read(_file)
    except Exception as _e1:
        if _verbose:
            print(f"{_e1}")
        raise Exception(f"Unable to read {_file}")
    else:
        if len(_candidates) <= 0:
            raise Exception(f"File {_file} contains no data")

    # connect to database
    try:
        engine = create_engine(DB_CONNECT)
        get_session = sessionmaker(bind=engine)
        session = get_session()
    except Exception as _e2:
        if _verbose:
            print(f"{_e2}")
        raise Exception(f"Failed to connect to database")

    _names = [_ for _ in _candidates['Name']]

    # gather co-ordinates into (Ra, Dec) tuples
    _coords = tuple(zip(_candidates['RA_deg'], _candidates['Dec_deg']))

    # Pass to milliquas_query
    qso, qprob, qoffset = milliquas_query(session, _coords, _names, _radius)
    print(qso,qprob,qoffset)
    session.close()

def milliquas_query(session, coords, names, _radius, _verbose: bool = True):

    # set variable(s)
    _begin = time.time()

    qso = []; qoffset = []; qprob = []

    # for all (RA, Dec) tuples, execute cone search and log candidates with qpct > 97(%)
    for _i, _e in enumerate(coords):

        # set up query
        try:
            # _radius /= 99999999999999
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
    print(f"Completed in {_end-_begin:.3f} sec")
    print(f"Found {len(qso)} QSOs in {len(coords)} candidates")

    return qso, qprob, qoffset

def asassn_query(session, coords, RADIUS_ARCSEC, _verbose: bool = False):

    _begin = time.time()

    f_asassn = 0.0 * len(coords)


# +
# main()
# -
if __name__ == '__main__':

    # get command line argument(s)
    _p = argparse.ArgumentParser(description=f'catalogs.py', formatter_class=argparse.RawTextHelpFormatter)
    _p.add_argument(f'--file', default=FILE, help="""File, default '%(default)s'""")
    _p.add_argument(f'--radius', default=RADIUS_ARCSEC, help="""Radius (arcsec), default %(default)s""")
    _p.add_argument(f'--verbose', default=False, action='store_true', help=f'if present, produce more verbose output')
    _a = _p.parse_args()

    # execute
    try:
        static_cats_query(_file=_a.file, _radius=float(_a.radius), _verbose=bool(_a.verbose))
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
