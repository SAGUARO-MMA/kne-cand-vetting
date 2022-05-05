from astropy.io import ascii
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sassy_src.models import *

# Import specific galaxy catalog filters on sassy_src

from typing import Optional
from astropy.coordinates import SkyCoord
from astropy import units as u

import argparse
import os
import time

DB_HOST = os.getenv('DB_HOST', 'localhost')
DB_NAME = os.getenv('DB_NAME', 'sassy')
DB_PASS = os.getenv('DB_PASS', None)
DB_PORT = os.getenv('DB_PORT', 5432)
DB_USER = os.getenv('DB_USER', 'sassy')

DB_CONNECT = f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
# FILE = os.path.abspath(os.path.expanduser("comb_master_FINAL.dat"))

# Chose radius to be equal to a 100 kpc offset at ~100 Mpc
RADIUS_ARCMIN = 3.183

def galaxy_search(RA: float, Dec: float, _radius: float = RADIUS_ARCMIN, _verbose: bool = False) -> Optional[list]:

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
