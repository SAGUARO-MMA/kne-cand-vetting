'''
Code to check that candidates are not minor planets

Most of this code has been modified from the SAGUARO pipeline for CSS
'''
import os, time
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord
import numpy as np
import ephem  # PyEphem module

def convert_mpcorb_to_monthly_catalog(filename_in, filename_out):
    raw_catalog_list = []
    f = open(filename_in, "r")
    for line in f:
        if line == "\n":
            continue  # in case there's a blank
        else:  # line in the original data file
            raw_catalog_list.append(line.rstrip())
    f.close()

    for n in range(100):
        if raw_catalog_list[n] == '-' * 160:
            start_line = n + 1
            # to define the start of the actual data table,
            # which comes after ~30 lines of header text

    cropped_catalog_list = []
    # crop off the header
    for n in range(len(raw_catalog_list) - start_line):
        cropped_catalog_list.append(raw_catalog_list[n + start_line])

    full_catalog = []

    for obj_mpc in cropped_catalog_list:
        abs_m_H = obj_mpc[8:14].strip()
        slope_G = obj_mpc[14:20].strip()
        epoch = obj_mpc[20:26].strip()
        mean_anomaly_M = obj_mpc[26:36].strip()
        peri = obj_mpc[37:47].strip()
        node = obj_mpc[48:58].strip()
        inclin = obj_mpc[59:69].strip()
        eccen = obj_mpc[70:80].strip()
        motion_n = obj_mpc[80:92].strip()
        a = obj_mpc[92:104].strip()
        unc_U = obj_mpc[105:107].strip()
        readable_designation = obj_mpc[166:194].strip()

        # MPC format has a "packed" date, allowing the epoch to be stored in
        # fewer digits. However, this must be converted to mm/dd/yyyy format
        # for XEphem.
        epoch_x = f'{int(epoch[3], 36):02d}/{int(epoch[4], 36):02d}.0/{int(epoch[0], 36):02d}{epoch[1:3]}'

        if unc_U == "":
            unc_U = "?"
        expanded_designation = readable_designation + " " + unc_U

        # Write XEphem format orbit to the full_catalog list.
        full_catalog.append(expanded_designation + ",e," + inclin + ","
                            + node + "," + peri + "," + a + "," + motion_n + "," + eccen + "," +
                            mean_anomaly_M + "," + epoch_x + "," + "2000" + ",H " + abs_m_H +
                            "," + slope_G + "\n")

    f2 = open(filename_out, "w")
    for obj in full_catalog:
        f2.write(obj)
    f2.close()

def movingobjectcatalog(obsmjd):
    catalog_list = []
    tobs = Time(obsmjd, format='mjd').strftime('%Y_%m')
    fnam = f"{tobs}_ORB.DAT"
    if not os.path.exists(fnam):
        os.system('wget -nv -O MPCORB.DAT http://www.minorplanetcenter.org/iau/MPCORB/MPCORB.DAT')
        convert_mpcorb_to_monthly_catalog('MPCORB.DAT', fnam)
    elif time.time() - os.path.getmtime(fnam)>(24*60*60):
        os.system('wget -nv -O MPCORB.DAT http://www.minorplanetcenter.org/iau/MPCORB/MPCORB.DAT')
        convert_mpcorb_to_monthly_catalog('MPCORB.DAT', fnam)
    with open(fnam) as f_catalog:
        for line in f_catalog:
            catalog_list.append(ephem.readdb(line))
    return catalog_list

def movingobjectfilter(s_catalog, s_ra, s_dec, obsmjd, filter_radius):
    """Searches for matches between (ra, dec, time) and the Minor Planet Center catalog. Takes about a minute to run."""
    tobs = Time(obsmjd, format='mjd')
    s_date = tobs.datetime
    ras, decs = [], []
    for body in s_catalog:
        body.compute(s_date)
        ras.append(body.a_ra)
        decs.append(body.a_dec)
    catalog_coords = SkyCoord(ras, decs, unit='rad')
    s_coords = SkyCoord(s_ra, s_dec, unit='deg')
    _, separation, _ = s_coords.match_to_catalog_sky(catalog_coords)
    return separation.arcsec < filter_radius

def is_minor_planet(ra, dec, discovery_mjd, filter_radius=25):
    """Checks if the target is a minor planet

    Args:
        ra (float|str) : float or string to pass to SkyCoord
        dec (float|str) : float or string declination to pass to SkyCoord
        discovery_mjd (float) : The discovery time to check for in MPC
        filter_radius (float) : Flag as a minor planet if less than this filter radius in arcseconds

    Returns:
        True if the candidate is within filter_radius at discovery_mjd, False otherwise
    """

    # download/load the moving object catalog from the Minor Planet Center
    print('Loading moving object catalog...')
    catalog = movingobjectcatalog(Time.now().mjd)
    print('Moving object catalog loaded.')
    
    return movingobjectfilter(catalog, ra, dec, discovery_mjd, filter_radius)
