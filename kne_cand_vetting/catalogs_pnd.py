#!/usr/bin/env python3


# +
# import(s)
# -
from astropy.io import ascii
from sassy_src.models.milliquas_q3c_orm import *
from sassy_src.models.milliquas_q3c_orm_filters import *
from typing import Optional

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
DB_CONNECT = f"postgresql+psycopg2://sassy:***REMOVED***@localhost:5432/sassy"
FILE = os.path.abspath(os.path.expanduser("comb_master_FINAL.dat"))
RADIUS_ARCSEC = 2.0


# +
# function: milliquas_query()
# -
#  ->  allows you to optionall accept inputs as lists
def milliquas_query(_file: str = FILE, _radius: float = RADIUS_ARCSEC, _verbose: bool = False) -> Optional[list]:
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

    # set variable(s)
    _begin = time.time()
    _names = [_ for _ in _candidates['Name']]
    qprob = [0.0] * len(_candidates)
    qso = []

    # gather co-ordinates into (Ra, Dec) tuples
    _coords = tuple(zip(_candidates['RA_deg'], _candidates['Dec_deg']))

    # for all (RA, Dec) tuples, execute cone search and log candidates with qpct > 97(%)
    for _i, _e in enumerate(_coords):

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
        for _x in MilliQuasQ3cRecord.serialize_list(query.all()):
            print(f'>>> QUASAR MATCH at RA, Dec = ({_e[0]},{_e[1]}), index={_i}!')

            # add the query dictionary to the qso list and modify the qprob list
            qso.append({**_x, **{'Candidate': _names[_i], 'Probability': 1.0, 'Candidate_RA': _e[0], 'Candidate_Dec': _e[1]}})
            qprob[_i] = 1.0

    # done
    session.close()
    _end = time.time()

    # return list of probabilities (although I think you'd be better off returning the qso list!)
    print(f"QSO Candidates: {qso}")
    print(f"QSO Probabilities: {qprob}")
    print(f"Completed in {_end-_begin:.3f} sec")
    print(f"Found {len(qso)} QSOs in {len(_candidates)} candidates")
    return qprob


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
        milliquas_query(_file=_a.file, _radius=float(_a.radius), _verbose=bool(_a.verbose))
    except Exception as _x:
        print(f"{_x}")
        print(f"Use:{__doc__}")
