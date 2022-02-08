import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

def milliquas_query(ra,dec,rad=2):

    '''
    Query the Million Quasar Catalog (Flesch 2021) for matches to
    kilonova candidates

    PARAMETERS
    ----------
    ra, dec : array of floats
        degrees
    radius : float
        default = 2, arcseconds

    RETURNS
    -------
    match : array of floats
        1. = match to milliquas source with qpct > 97; 0. = no match

    '''
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    print(coords)

    # Query statement --> recall you're using an array

    # Only accept quasars with Qpct > 97.
    # qpct = r['Qpct']
    q =  np.zeros(len(ra))
    # q[qpct > 97.] = 1.0

    return coords # q
