# JCR May 2022, adapted from Charlie Kilpatrick, Vic Dong et al.

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import numpy as np

from sassy_q3c_models.glade_plus_q3c_orm import GladePlusQ3cRecord
from sassy_q3c_models.glade_plus_q3c_orm_filters import glade_plus_q3c_orm_filters
from sassy_q3c_models.gwgc_q3c_orm import GwgcQ3cRecord
from sassy_q3c_models.gwgc_q3c_orm_filters import gwgc_q3c_orm_filters
from sassy_q3c_models.hecate_q3c_orm import HecateQ3cRecord
from sassy_q3c_models.hecate_q3c_orm_filters import hecate_q3c_orm_filters
from sassy_q3c_models.sdss12photoz_q3c_orm import Sdss12PhotoZQ3cRecord
from sassy_q3c_models.sdss12photoz_q3c_orm_filters import sdss12photoz_q3c_orm_filters
from sassy_q3c_models.ps1_q3c_orm import Ps1Q3cRecord
from sassy_q3c_models.ps1_q3c_orm_filters import ps1_q3c_orm_filters
from sassy_q3c_models.desi_spec_q3c_orm import DesiSpecQ3cRecord
from sassy_q3c_models.desi_spec_q3c_orm_filters import desi_spec_q3c_orm_filters
from sassy_q3c_models.ls_dr10_q3c_orm import LsDr10Q3cRecord
from sassy_q3c_models.ls_dr10_q3c_orm_filters import ls_dr10_q3c_orm_filters

from typing import Optional
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=69.6 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
c_over_H0 = (const.c / cosmo.H0).to_value(u.Mpc)

import argparse
import os
import time

from urllib.request import urlopen
import json

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
RADIUS_ARCMIN = 3.437745624870037
PCC_THRESHOLD = 0.80

def galaxy_search(RA: float, Dec: float, _radius: float = RADIUS_ARCMIN, _pcc_thresh: float = PCC_THRESHOLD,
                  _verbose: bool = False, db_connect: str = DB_CONNECT) -> Optional[list]:

    """
    Searches for galaxy match to candidate in GLADE, HECATE, GWGC, SDSS, Legacy DR8 and PS1-STRM catalogs
    Returns 10 galaxies with lowest Pcc values (must be < PCC_THRESHOLD), their offsets, redshifts and magnitudes

    PARAMETERS
    ----------
    RA, Dec : transient coordinates, float
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
    dist : spectroscopic redshift, -99.0 if none, list of floats
    z : photometric redshift, -99.0 if none, list of floats
    etc.
    """

    _radius /= 60.0 # convert to degrees

    # connect to database
    try:
        engine = create_engine(db_connect)
        get_session = sessionmaker(bind=engine)
        session = get_session()
    except Exception as _e2:
        if _verbose:
            print(f"{_e2}")
        raise Exception(f"Failed to connect to database")

    # loop through RA, Dec here
    _begin = time.time()
    glade=0; gwgc=0; hecate=0; sdss=0; ps1=0; desi=0; lsdr10=0

    # Find matches in GLADE:
    GLADE_matches, GLADE_ra, GLADE_dec, GLADE_offset, GLADE_mag, GLADE_filt, GLADE_dist, GLADE_dist_err, GLADE_distflag, GLADE_source, GLADE_name, GLADE_z, GLADE_zerr = query_GLADE(session, RA, Dec, _radius)
    if GLADE_matches>0:
        glade+=1
        
    # Find matches in GWGC:
    GWGC_matches, GWGC_ra, GWGC_dec, GWGC_offset, GWGC_mag, GWGC_filt, GWGC_dist, GWGC_dist_err, GWGC_source, GWGC_name = query_GWGC(session, RA, Dec, _radius)

    if GWGC_matches>0:
        gwgc+=1    
    
    # convert luminosity distance to redshift
    GWGC_z = np.array(GWGC_dist) / c_over_H0
    GWGC_zerr = np.array(GWGC_dist_err) / c_over_H0
        
    HECATE_matches, HECATE_ra, HECATE_dec, HECATE_offset, HECATE_mag, HECATE_filt, HECATE_dist, HECATE_dist_err, HECATE_source, HECATE_name, HECATE_v, HECATE_verr = query_hecate(session, RA, Dec, _radius)

    if HECATE_matches>0:
        hecate+=1
        
    # convert recession velocity to redshift
    HECATE_z = np.array(HECATE_v) / const.c.to_value(u.km / u.s)
    HECATE_zerr = np.array(HECATE_verr) / const.c.to_value(u.km / u.s)
        
    DESI_matches, DESI_ra, DESI_dec, DESI_offset, DESI_mag, DESI_filt, DESI_z, DESI_zerr, DESI_source, DESI_name = query_desi_spec(session, RA, Dec, _radius)

    if DESI_matches>0:
        desi+=1

    # convert redshift to distance
    DESI_dist = np.array(DESI_z) * c_over_H0
    DESI_dist_err = np.array(DESI_zerr) * c_over_H0
        
    SDSS_matches, SDSS_ra, SDSS_dec, SDSS_offset, SDSS_mag, SDSS_filt, SDSS_z, SDSS_zerr, SDSS_source, SDSS_name = query_sdss12phot(session, RA, Dec, _radius)

    if SDSS_matches>0:
        sdss+=1

    # convert redshift to distance
    SDSS_dist = np.array(SDSS_z) * c_over_H0
    SDSS_dist_err = np.array(SDSS_zerr) * c_over_H0
        
    PS1_matches, PS1_ra, PS1_dec, PS1_offset, PS1_mag, PS1_filt, PS1_z, PS1_zerr, PS1_source, PS1_name = query_ps1(session, RA, Dec, _radius)

    if PS1_matches>0:
        ps1+=1

    # convert redshift to distance
    PS1_dist = np.array(PS1_z) * c_over_H0
    PS1_dist_err = np.array(PS1_zerr) * c_over_H0
    
    LSDR10_matches, LSDR10_ra, LSDR10_dec, LSDR10_offset, LSDR10_mag,LSDR10_filt, LSDR10_z, LSDR10_zerr, LSDR10_source, LSDR10_name = query_LS_DR10_photoz(session, RA, Dec, _radius)

    if LSDR10_matches>0:
        lsdr10+=1

    LSDR10_dist = np.array(LSDR10_z) * c_over_H0
    LSDR10_dist_err = np.array(LSDR10_zerr, dtype='object') * c_over_H0
    LSDR10_zerr = np.array(LSDR10_zerr, dtype='object')

    # sum the findings, turn into numpy arrays
    tot_names = np.array(GLADE_name + GWGC_name + HECATE_name + DESI_name + SDSS_name + PS1_name + LSDR10_name, dtype=str)
    tot_offsets = np.array(GLADE_offset + GWGC_offset + HECATE_offset + DESI_offset + SDSS_offset + PS1_offset + LSDR10_offset)
    tot_mags = np.array(GLADE_mag + GWGC_mag + HECATE_mag + DESI_mag + SDSS_mag + PS1_mag + LSDR10_mag)
    tot_ra = np.array(GLADE_ra + GWGC_ra + HECATE_ra + DESI_ra + SDSS_ra + PS1_ra + LSDR10_ra)
    tot_dec = np.array(GLADE_dec + GWGC_dec + HECATE_dec + DESI_dec + SDSS_dec + PS1_dec + LSDR10_dec)
    tot_filt = np.array(GLADE_filt + GWGC_filt + HECATE_filt + DESI_filt + SDSS_filt + PS1_filt + LSDR10_filt)
    tot_dists = np.concatenate([GLADE_dist, GWGC_dist, HECATE_dist, DESI_dist, SDSS_dist, PS1_dist, LSDR10_dist])
    tot_z = np.concatenate([GLADE_z, GWGC_z, HECATE_z, DESI_z, SDSS_z, PS1_z, LSDR10_z])
    tot_source = np.array(GLADE_source + GWGC_source + HECATE_source + DESI_source + SDSS_source + PS1_source + LSDR10_source)

    # also sum the error findings, they require some special treatment though
    tot_dist_errs = np.array(
        list(GLADE_dist_err) +
        list(GWGC_dist_err) +
        list(HECATE_dist_err) +
        list(DESI_dist_err) +
        list(SDSS_dist_err) +
        list(PS1_dist_err) +
        list(LSDR10_dist_err),
        dtype='object'
    )

    tot_zerr = np.array(
        list(GLADE_zerr) +
        list(GWGC_zerr) +
        list(HECATE_zerr) +
        list(DESI_zerr) +
        list(SDSS_zerr) +
        list(PS1_zerr) +
        list(LSDR10_zerr),
        dtype='object'
    )

    
    PCCS = pcc(tot_offsets,tot_mags)
    
    # put some basic cut on Pcc ?
    pcc_args = np.argsort(PCCS)[:10]
    cond = (PCCS[pcc_args] < _pcc_thresh)

    _end = time.time()

    print(f"Completed galaxy search in {_end-_begin:.3f} sec")
    if glade==1:
        print(f"Found a GLADE galaxy match.")
    if gwgc==1:
        print(f"Found a GWGC galaxy match.")
    if hecate==1:
        print(f"Found a HECATE galaxy match.")
    if sdss==1:
        print(f"Found SDSS DR12 Photo-z Catalog galaxy match.")
    if lsdr10==1:
        print(f"Found Legacy Survey DR10 Photo-z Catalog galaxy match.")

    all_data = [
        {
            'ID':tot_names[pcc_args][cond][i],
            'PCC':PCCS[pcc_args][cond][i],
            'Offset':tot_offsets[pcc_args][cond][i],
            'RA':tot_ra[pcc_args][cond][i],
            'Dec':tot_dec[pcc_args][cond][i],
            'Dist':tot_dists[pcc_args][cond][i],
            # need to convert all to a np array first so we can convert to a list
            # this is cause np arrays are not json serializable
            'DistErr':np.array(tot_dist_errs[pcc_args][cond][i]).tolist(), 
            'z':tot_z[pcc_args][cond][i],
            # same thing is true here for zErr
            'zErr':np.array(tot_zerr[pcc_args][cond][i]).tolist(),
            'Mags':tot_mags[pcc_args][cond][i],
            'Filter':tot_filt[pcc_args][cond][i],
            'Source':tot_source[pcc_args][cond][i]
        } for i in range(len(PCCS[pcc_args][cond]))
    ]

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
    sigma = (1/(0.33*np.log(10)))*10**(0.33*(m-24)-2.44)
    prob = 1-np.exp(-(np.pi*(r**2)*sigma))

    return prob

def nanomgy_to_mag(nmgy):
    mag = 22.5 - 2.5 * np.log10(nmgy)
    return mag

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
            elif key=='hyperleda':
                if key.startswith('SDSS')==False:
                    name = key+_dict[key]
            else:
                print(key,_dict[key])
                name = _dict[key]
            break
        if key==keys[-1]:
            name = catalog+str(_dict[key])

    return name

def query_LS_DR10_photoz(session, ra, dec, _radius, _verbose: bool = True):
    """
    Query Legacy Survey DR 10 Photo-z Catalog
    """
    m=0
    gal_offset = []; mag = []; filt = []; z = []; z_err = []; gal_ra = [];
    gal_dec = []; source = []; name = []

    try:
        query = session.query(LsDr10Q3cRecord)
        query = ls_dr10_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}'})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    if len(query.all()) > 0:
        m+=1
        for _x in LsDr10Q3cRecord.serialize_list(query.all()):
            if np.isfinite(_x['flux_r']) and _x['flux_r'] != -99:
                if _x['z_spec'] != -99:
                    z.append(_x['z_spec'])
                    z_err.append(0.)  # no error for spectroscopic redshift
                elif _x['z_phot_mean'] != -99:
                    z.append(_x['z_phot_mean'])

                    # tuple of lower and upper errorbars cause LS_DR10 has both instead
                    # of a single error. We will catch this later.
                    z_err.append(np.array([_x['z_phot_l68'],_x['z_phot_u68']])) 
                else:
                    continue
                mag.append(nanomgy_to_mag(_x['flux_r']))
                filt.append('r')
                gal = SkyCoord(_x['ra']*u.deg, _x['declination']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['ra'])
                gal_dec.append(_x['declination'])
                source.append('LS_DR10')
                name.append(_x['lid'])

    return m, gal_ra, gal_dec, gal_offset, mag, filt, z, z_err, source, name

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
    gal_offset = []; mag = []; filt = []; dist = []; dist_err = []; gal_ra = []; gal_dec = []; distflag = []; source = []; name = []
    z = []
    z_err = []

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
                z.append(_x['z_helio'])
                z_err.append(_x['z_err'])
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
                z.append(_x['z_helio'])
                z_err.append(_x['z_err'])

    return m, gal_ra, gal_dec, gal_offset, mag, filt, dist, dist_err, distflag, source, name, z, z_err

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
    v = []
    v_err = []

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
                v.append(_x['v'])
                v_err.append(_x['e_v'])

    return m, gal_ra, gal_dec, gal_offset, mag, filt, dist, dist_err, source, name, v, v_err

def query_desi_spec(session, ra, dec, _radius, _verbose: bool = True):

    # Query the DESI spectroscopic catalog
    m=0
    gal_offset = []; mag = []; filt = []; z = []; z_err = []; gal_ra = [];
    gal_dec = []; source = []; name = []

    # set up query
    try:
        query = session.query(DesiSpecQ3cRecord)
        query = desi_spec_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}'})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    if len(query.all()) > 0:
        m+=1
        for _x in DesiSpecQ3cRecord.serialize_list(query.all()):
            if _x['z']>0:
                if _x['flux_r']>0:
                    mag_r = - 2.5 * np.log10(_x['flux_r']*10**-9) # convert nmy to Jy
                    mag.append(mag_r)
                    filt.append('r')
                elif _x['gaia_phot_g_mean_mag']>0:
                    mag.append(_x['gaia_phot_g_mean_mag'])
                    filt.append('G')
                else:
                    continue
                z.append(_x['z'])
                z_err.append(_x['zerr'])
                gal = SkyCoord(_x['target_ra']*u.deg, _x['target_dec']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['target_ra'])
                gal_dec.append(_x['target_dec'])
                source.append('DESI')
                name.append(_x['targetid'])

    return m, gal_ra, gal_dec, gal_offset, mag, filt, z, z_err, source, name


def query_sdss12phot(session, ra, dec, _radius, _verbose: bool = True):

    # Query the SDSS DR12 Photo-z Catalog

    m=0
    gal_offset = []; mag = []; filt = []; z = []; z_err = []; gal_ra = [];
    gal_dec = []; source = []; name = []

    try:
        query = session.query(Sdss12PhotoZQ3cRecord)
        query = sdss12photoz_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}'})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    if len(query.all()) > 0:
        m+=1
        for _x in Sdss12PhotoZQ3cRecord.serialize_list(query.all()):
            if np.isfinite(_x['rmag']) and _x['rmag'] != -9999.:
                if np.isfinite(_x['zsp']) and _x['zsp'] != -9999.:
                    z.append(_x['zsp'])
                    z_err.append(0.)  # no error for spectroscopic redshift
                elif np.isfinite(_x['zph']) and _x['zph'] != -9999.:
                    z.append(_x['zph'])
                    z_err.append(_x['e_zph'])
                else:
                    continue
                mag.append(_x['rmag'])
                filt.append('r')
                gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['ra'])
                gal_dec.append(_x['dec'])
                source.append('SDSS_DR12')
                name.append(_x['sdss12'])

    return m, gal_ra, gal_dec, gal_offset, mag, filt, z, z_err, source, name

def query_ps1(session, ra, dec, _radius, _verbose: bool = True):

    # Query the PS1 STRM Photo-z Catalog

    gal_offset = []; mag = []; filt = []; z = []; z_err = []; gal_ra = [];
    gal_dec = []; source = []; name = []; m=0

    try:
        query = session.query(Ps1Q3cRecord)
        query = ps1_q3c_orm_filters(query, {'cone': f'{ra},{dec},{_radius}', 'ps_score__lte': 0.83})
    except Exception as _e3:
        if _verbose:
            print(f"{_e3}")
        print(f"Failed to execute query for RA, Dec = ({ra},{dec})")

    if len(query.all()) > 0:
        m+=1
        for _x in Ps1Q3cRecord.serialize_list(query.all()):

            #### DO NOT HAVE MAGNITUDE YET
            if _x['rmeanpsfmag'] is not None and _x['rmeanpsfmag'] != -999.:
                if np.isfinite(_x['z_phot']) and _x['z_phot'] != -999.:
                    z.append(_x['z_phot'])
                    z_err.append(_x['z_err'])
                else:
                    continue
                mag.append(_x['rmeanpsfmag'])
                filt.append('r')
                gal = SkyCoord(_x['ra']*u.deg, _x['dec']*u.deg)
                cand = SkyCoord(ra*u.deg, dec*u.deg)
                gal_offset.append(cand.separation(gal).arcsec)
                gal_ra.append(_x['ra'])
                gal_dec.append(_x['dec'])
                source.append('PS1_STRM')
                name.append(_x['psps_objid'])

    return m, gal_ra, gal_dec, gal_offset, mag, filt, z, z_err, source, name
    
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
        # engine = create_engine(DB_CONNECT)
        # get_session = sessionmaker(bind=engine)
        # session = get_session()
        output = galaxy_search(RA=_a.RA, Dec =_a.Dec, _radius=float(_a.radius), _pcc_thresh=float(_a.pcc_threshold), _verbose=bool(_a.verbose))
        # output = query_sdss12phot(session, _a.RA[0], _a.Dec[0], _radius=float(_a.radius), _verbose=bool(_a.verbose))
        print('Finished.')
        print(output)
    except Exception as _x:
        print(f"{_x}")
        print(f"Use:{__doc__}")
