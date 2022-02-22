import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, Float, union, text, and_
from sqlalchemy.sql import select
from astropy.io import fits, ascii
import time, os, sys

import sqlalchemy
print(sqlalchemy.__version__)

# ssh sassy2
# PGPASSWORD=***REMOVED*** psql --echo-all -h localhost -p 5432 -U sassy -d sassy -c "SELECT * FROM milliquas_q3c;"

def milliquas_query(file,rad=2):

    '''
    Query the Million Quasar Catalog (Flesch 2021) for matches to
    kilonova candidates

    PARAMETERS
    ----------
    ra, dec : float --> right now it's an ascii file that I pull from
        degrees
    radius : float
        default = 2, arcseconds

    RETURNS
    -------
    match : float
        1. = match to milliquas source with qpct > 97; 0. = no match

    '''
    # set up parameters to connect to sassy database where milliquas_q3c table stored
    dbname = 'sassy'
    eng = 'postgresql://sassy:***REMOVED***@localhost:5432/' + dbname

    engine = create_engine(eng, echo = True)
    meta = MetaData()

    conn = engine.connect()

    # set up what the columns will look like
    milliquas = Table(
    'milliquas_q3c', meta,
    Column('ra', Float),
    Column('dec', Float),
    Column('qpct',Float)
    )

    all_cand = ascii.read(file)
    tot_names = all_cand['Name']
    tot_radeg = all_cand['RA_deg'].astype(float)
    tot_decdeg = all_cand['Dec_deg'].astype(float)

    t0 = time.time()
    qso = []
    for i in range(len(tot_radeg)):
        ra = tot_radeg[i]
        dec = tot_decdeg[i]

        # need to make conditions but now this should work!!
        ra_high = ra + (rad / 3600.)
        ra_low = ra - (rad / 3600.)
        dec_high = dec + (rad / 3600.)
        dec_low = dec - (rad / 3600.)

        # set conditions to select match
        s = milliquas.select().where(and_(text("milliquas_q3c.ra between :x1 and :x2"),text("milliquas_q3c.dec between :y1 and :y2"),text("milliquas_q3c.qpct>97.")))

        result = conn.execute(s, x1 = ra_low, x2 = ra_high, y1 = dec_low, y2 = dec_high).fetchall()

        if len(result) > 0:
            print('QUASAR MATCH!')
            qso.append(1.0)
        else:
            print('No quasar match.')
            qso.append(0.0)

    print("Completed in {:.3f} sec".format(time.time()-t0))
    print("Found {:.1f} QSO of {:.2f} candidates".format(sum(qso),len(qso)))

milliquas_query('comb_master_FINAL.dat')

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
