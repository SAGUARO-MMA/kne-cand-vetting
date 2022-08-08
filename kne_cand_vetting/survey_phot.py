### Query public surveys for photometric detections prior to GW event
### Currently working on ZTF only

from astropy.coordinates import SkyCoord
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from astropy.time import Time
import os
import requests
import argparse
from datetime import datetime, timedelta

from sassy_q3c_models.ztf_q3c_orm import ZtfQ3cRecord
from sassy_q3c_models.ztf_q3c_orm_filters import *
from sassy_q3c_models.ztf_q3c_orm_cli import *

# constants -- are these still needed?
DB_HOST = os.getenv('POSTGRES_HOST', 'localhost')
DB_NAME = os.getenv('POSTGRES_DB', 'sassy')
DB_PASS = os.getenv('POSTGRES_PASSWORD', None)
DB_PORT = os.getenv('POSTGRES_PORT', 5432)
DB_USER = os.getenv('POSTGRES_USER', 'sassy')

DB_CONNECT = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

RADIUS_ARCSEC = 2.0

def query_ZTFpubphot(RA: float, Dec: float, _radius: float = RADIUS_ARCSEC, _verbose: bool = False):

    _radius /= 3600 # converting to Degrees
    if _radius <= 0.0:
        raise Exception(f"Invalid input, _radius={_radius}")

    try:
        engine = create_engine(DB_CONNECT)
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

    short_keys = ['jd', 'magpsf', 'sigmapsf', 'filtername', 'diffmaglim']
    ztfdict = ZtfQ3cRecord.serialize_list(query.all())
    short_ztfdict = {key: [det['candidate'][key] for det in ztfdict] for key in short_keys}

    if len(short_ztfdict)>0:
        print('{0} photometric detections found in ZTF.'.format(len(short_ztfdict)))

    return short_ztfdict

def ATLAS_forcedphot(t_Event: datetime, RA: float, Dec: float, _radius: float = RADIUS_ARCSEC, _verbose: bool = False):

    t_querystart = t_Event - timedelta(days=200)
    t_queryend = datetime.now()


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
        # For ZTF query entire database
        ztfphot = query_ZTFpubphot(RA=a.RA, Dec =a.Dec, _radius=float(a.radius), _verbose=bool(a.verbose))
        # For ATLAS + other forced phot, only query -200 days
    except Exception as _x:
        print(f"{_x}")
        print(f"Use:{__doc__}")
