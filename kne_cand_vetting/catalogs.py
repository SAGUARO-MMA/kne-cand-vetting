import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, union
import pandas, time

# ssh sassy2
# PGPASSWORD=***REMOVED*** psql --echo-all -h localhost -p 5432 -U sassy -d sassy -c "SELECT * FROM milliquas_q3c;"
# -d: database
# -h: host
# -U: user (uid)

def SQL_query(dbname):

    # filled these in correctly?
    # dialect[+driver]://+ dsn_uid + ':' + dsn_pwd + '@'+dsn_hostname+':'+dsn_port+'/' + dsn_database
    'PostgreSQL://+sassy:***REMOVED***@localhost:5432/' + dbname

    engine = create_engine(*args, echo = True)
    meta = MetaData()

    conn = engine.connect()
    sql = "SELECT * FROM table1"
    table1 = # get
    u = union(addresses.select().where(addresses.c.email_add.like('%@gmail.com addresses.select().where(addresses.c.email_add.like('%@yahoo.com'))))

    t0 = time.time()
    df = pd.read_sql_query(sql, engine)
    print("Completed in {:.1f} sec".format(time.time()-t0))

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

    query = """
    select p.*
    inner join fGetNearbyObjEq({ra},{dec},{radius}) as r on p.objid=r.objid
    """.format(**locals())

    # Only accept quasars with Qpct > 97.
    # qpct = r['Qpct']
    q =  np.zeros(len(ra))
    # q[qpct > 97.] = 1.0

    return coords # q

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
