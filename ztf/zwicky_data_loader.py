import argparse
import io
import os
import re
import astroquery.exceptions
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astroquery.ipac.ned import Ned
import logging

from config import *

logging.basicConfig(level=logging.INFO)

def format_to_3e(x):
    """
    Format a number to 3 decimal scientific notation.

    Parameters
    ----------
    x : float
        The number to format.

    Returns
    -------
    str
        The number formatted in scientific notation with 3 decimal places.
    """
    return "{:.3e}".format(x)


def query_sky_coord(source_name):
    """
    Query the SkyCoord of a source.

    Parameters
    ----------
    source_name : str
        The name of the source.

    Returns
    -------
    SkyCoord
        The SkyCoord object for the source.
    """
    coords = SkyCoord.from_name(source_name)
    return coords


def calculate_norm(flux, flux_err):
    """
    Calculate the normalization factor for flux and flux error.

    Parameters
    ----------
    flux : pd.Series
        The flux.
    flux_err : pd.Series
        The flux error.

    Returns
    -------
    int, pd.Series, pd.Series
        The log of the normalization factor, the normalized flux, and the normalized flux error.
    """
    if flux.size > 0:
        avr = flux.mean()
        log = int(math.floor(math.log10(avr)))
        norm = 10 ** log
        flux = flux / norm
        flux_err = flux_err / norm
        return log, flux, flux_err
    else:
        return 0, flux, flux_err


def query_simbad(source_name):
    """
    Query the coordinates of a source using Simbad.

    Parameters
    ----------
    source_name : str
        The name of the source.

    Returns
    -------
    SkyCoord
        The SkyCoord object for the source.
    """
    from astroquery.simbad import Simbad
    simbad_queried_object = Simbad.query_object(source_name)
    if not simbad_queried_object:
        return query_sky_coord(source_name)
    else:
        coords = SkyCoord(f"{simbad_queried_object['RA'][0]} {simbad_queried_object['DEC'][0]}",
                          unit=("hourangle", "deg"))
        return coords


def query_coordinates(source_name):
    """
    This function queries the source name to NED, SIMBAD or SkyCoord databases and returns coordinates of that source.

    Parameters
    ----------
    source_name :   str
        Name of the source for what you want to find coordinates.

    Returns
    -------
    astropy.coordinates.sky_coordinate.SkyCoord
        Coordinates of the source.
    """
    try:
        ned_queried_object = Ned.query_object(source_name)
        if ned_queried_object:
            coords = SkyCoord(f"{ned_queried_object['RA'][0]} {ned_queried_object['DEC'][0]}", unit="deg")
        else:
            coords = query_simbad(source_name)
    except astroquery.exceptions.RemoteServiceError:
        coords = query_simbad(source_name)

    if coords:
        res_coords = SkyCoord(f"{math.round(coords.ra.value, 4)} {math.round(coords.dec.value, 4)}", unit="deg")  # noqa
        return res_coords
    else:
        print("NED, SIMBAD and SkyCoord queries failed for the requested source name. \n "
              "Please try rerun with coordinates!")
        raise Exception("NED, SIMBAD and SkyCoord queries failed for the requested source name. \n"
                        "Please try rerun with coordinates! ")


def query_ztf_data(coordinates):
    """
    Query ZTF light curve data for a given set of sky coordinates.

    Parameters
    ----------
    coordinates : `astropy.coordinates.SkyCoord`
        The sky coordinates for which ZTF data is queried.

    Returns
    -------
    data_df : `pandas.DataFrame`
        A `pandas` data frame containing the queried ZTF light curve data.

    Raises
    ------
    Exception
        If there are no Zwicky observations for the source, an exception is raised with the message
        "There are no Zwicky observations for this source."
    """
    prefix_of_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"
    url = f"{prefix_of_url}POS=CIRCLE%20{coordinates.ra.deg}%20{coordinates.dec.deg}%200.0014&BANDNAME=g,r,i&FORMAT=csv"
    requested_data = requests.get(url).content
    # noinspection PyTypeChecker
    data_df = pd.read_csv(io.StringIO(requested_data.decode('utf-8')))
    if data_df.size == 0:
        print('There are no Zwicky observations for this source.')
        raise Exception('There are no Zwicky observations for this source.')
    return data_df, url


def calc_nh_value_requests(coordinates):
    """
    Calculate the NH value for a given set of sky coordinates.

    Parameters
    ----------
    coordinates : `astropy.coordinates.SkyCoord`
        The sky coordinates for which nh value is to be calculated.

    Returns
    -------
    nh_value : float
        The calculated nh value.
    """
    entry = f"{coordinates.ra.value},{coordinates.dec.value}"

    url = "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl"

    params = {
        'Entry': entry,  # Replace with your object name or coordinates
        'NR': 'GRB/SIMBAD+Sesame/NED',
        'CoordSys': 'Equatorial',
        'equinox': '2000',
        'radius': '0.1',
        'usemap': '0'
    }

    # Make a GET request with the specified parameters
    response = requests.get(url, params=params)
    # Check the response status code
    if response.status_code == 200:
        pattern = r"Weighted average nH \(cm\*\*-2\)\s+(\d+\.\d+E[+-]\d+)"
        match = re.search(pattern, response.text)

        if match:
            nh_value = match.group(1)
            return float(nh_value)
        else:
            print("nH value not found in the response.")
            return None
    else:
        print(f"Request failed with status code: {response.status_code}")
        return None


def calc_nh_value(coordinates):
    """
    Calculate the NH value for a given set of sky coordinates.

    Parameters
    ----------
    coordinates : `astropy.coordinates.SkyCoord`
        The sky coordinates for which nh value is to be calculated.

    Returns
    -------
    nh_value : float
        The calculated nh value.
    """
    os.system(f"nh equinox=2000 ra={coordinates.ra.deg} dec={coordinates.dec.deg} > nh_value_tmp.txt")
    nh_value = os.popen("echo `tail -1 nh_value_tmp.txt |  awk '{print $7}'`").read().strip()
    os.system('rm nh_value_tmp.txt')
    return float(nh_value)


def find_lambda_and_const(filter_name):
    """
    Find lambda and constant values for a given filter.

    Parameters
    ----------
    filter_name : str
        The name of the filter for which lambda and constant values are to be found.

    Returns
    -------
    lambda_ : float
        The lambda value for the given filter.
    const_ : float
        The constant value for the given filter.

    Raises
    ------
    ValueError
        If the filter is not supported, a ValueError is raised with the message "mag2flux: Filter not supported".
    """
    if filter_name not in list_of_filters:
        print("mag2flux: Filter not supported")
        raise ValueError("mag2flux: Filter not supported")
    else:
        lambda_, const_ = filter_data.get(filter_name, (None, None))
    if not lambda_ and not const_:
        print('mag2flux: Filter not supported')
        return
    else:
        return lambda_, const_


def calculate_a_band(lambda_, nh_value):
    """
    Calculates the A-band value based on the input wavelength and NH value.

    Parameters
    ----------
    lambda_ : float
        Wavelength in Angstroms.
    nh_value : float
        NH value in cm^-2.

    Returns
    -------
    a_band: float
        The calculated A-band value.
    """
    x = (10000 / lambda_)
    av = max(Rv * (-0.055 + nh_value * 1.987e-22), 0)
    ebv = av / Rv
    if 1.1 > x >= 0.3:
        a_band = (0.574 * (x ** 1.61) - 0.527 * (x ** 1.61) / Rv) * av
    elif 3.3 > x >= 1.1:
        y = x - 1.82
        aa = 1 + (0.17699 * y) - (0.50447 * y ** 2) - (0.02427 * y ** 3) + (0.72085 * y ** 4) + (0.01979 * y ** 5) - (
                0.77530 * y ** 6) + (0.32999 * y ** 7)
        bb = 1.41338 * y + (2.28305 * y ** 2) + (1.07233 * y ** 3) - (5.38434 * y ** 4) - (0.62251 * y ** 5) + (
                5.30260 * y ** 6) - (2.09002 * y ** 7)
        a_band = (aa + (bb / Rv)) * av
    elif 10.0 >= x >= 3.3:
        c2 = -0.824 + (4.717 / Rv)
        c1 = 2.03 - (3.007 * c2)
        dx = (x * x) / (((x ** 2) - (4.596 ** 2)) ** 2 + ((x * 0.99) ** 2))
        px = (0.5392 * ((x - 5.9) ** 2)) + (0.05644 * ((x - 5.9) ** 3))
        if x <= 5.9:
            px = 0
        a_band = (c1 + c2 * x + 3.23 * dx + 0.41 * px) * ebv + av
    else:
        a_band = 0
    return a_band


def calculate_flux_and_frequency(nh, filter_name, magnitude):
    """
    Calculates the flux and frequency based on the input filter name, magnitude, and NH value.

    Parameters
    -------
    nh : float
        NH value in cm^-2.
    filter_name : str
        Name of the filter.
    magnitude : float
        Magnitude.

    Returns
    -------
    pandas.Series: The calculated flux and frequency in a Series object.
    """
    lambda_, const = find_lambda_and_const(filter_name)
    a_band = calculate_a_band(lambda_, nh)
    frequency = c / (lambda_ * 1.e-8)
    if magnitude == 0:
        flux = 0
    else:
        flux = 10. ** (-0.4 * (magnitude - a_band) + const) * frequency
    return pd.Series([flux, frequency])


def process_results_df(data_in_df, output_file, url):
    logging.info("Starting process_results_df function...")

    try:
        # Replacing applymap with apply and lambda-map combination
        data_in_df[['frequency', 'flux', 'flux_err']] = data_in_df[['frequency', 'flux', 'flux_err']].apply(
            lambda col: col.map(format_to_3e))

        logging.info("Data formatting completed.")

        data_in_df['mjd_start'] = data_in_df['mjd'].round(6)
        data_in_df['Flag'] = ''
        data_in_df['Catalog'] = "ZTF"
        data_in_df['Reference'] = "Bellm et al. 2019 PASP 131 018002"
        data_in_df.columns = ["freq. ", "flux ", "err_flux ", "MJD_start ", "MJD_end ", "flag", "catalog", "reference"]
        data_in_df.to_csv(output_file, index=False)

        logging.info(f"Data processed and saved to: {output_file}")

    except TypeError as te:
        logging.error(f"TypeError encountered: {te}")
        logging.error(
            "This typically happens due to incorrect indexing. Check the function format_to_3e or dataframe structures.")

    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")


def run_with_coordinates(coords, path, output_file_name=None):
    """
    This function queries ZTF data for a given set of coordinates,
    calculates flux and frequency values based on the given coordinates,
    plots the results and saves the data to an output file.

    Parameters
    ----------
    coords : `astropy.coordinates.SkyCoord`
        Object containing the RA and Dec of the target source.
    path : `str`
        Path to directory where the output file will be saved.
    output_file_name : `str`, optional
        Name of the output file, by default None.

    Raises
    ------
    Exception
        If the given coordinates have Dec < -31.
    """
    if coords.dec.value < -31:
        print("The position is not covered by ZTF in case Dec < -31. ")
        raise Exception("The position is not covered by ZTF in case Dec < -31. ")
    data_df, url = query_ztf_data(coords)
    nh_value = calc_nh_value_requests(coords)
    if not output_file_name:
        output_file = "/VOU_Blazars/ZTF_data.csv"
    else:
        output_file = f"{path}/{output_file_name.replace(' ', '').replace('+', 'P').replace('-', 'M')}"
    data_df['filtercode'] = "ps" + data_df['filtercode'].str[1]
    data_df = pd.concat([data_df,
                         data_df.apply(lambda x: calculate_flux_and_frequency(nh_value, x['filtercode'], x['mag']),
                                       axis=1).rename(columns={0: 'flux', 1: 'frequency'})], axis=1)
    data_df = pd.concat([data_df,
                         data_df.apply(
                             lambda x: calculate_flux_and_frequency(nh_value, x['filtercode'], x['mag'] + x['magerr']),
                             axis=1).rename(columns={0: 'flux_err_minus', 1: 'frequency_minus'})], axis=1)
    data_df['flux_err'] = data_df.flux - data_df.flux_err_minus
    data = data_df.loc[:, ['frequency', 'flux', 'flux_err', 'mjd']]
    process_results_df(data, output_file, url)


def run_with_source_name(src_name, path):
    """
    This function queries the coordinates of a source based on its name.

    Parameters
    ----------
    src_name : str
        Name of the target source.
    path : str
        Path to directory where the output file will be saved.

    """
    coords = query_coordinates(src_name)
    run_with_coordinates(coords, path, src_name)


def run_for_sources(list_of_source_names, path):
    """
    Run pipeline for a list of source names.

    Parameters
    ----------
    list_of_source_names : list
        List of source names.
    path : str
        Path to save the output.

    """

    for src_name in list_of_source_names:
        try:
            run_with_source_name(src_name, path)
        except Exception:
            print(f'Data could not been downloaded for the source {src_name}.')
            continue
        else:
            print(f'Data for the source {src_name} have been downloaded.')


def read_source_names_from_file(source_names_filename, path):
    """
    Read source names from file and run pipeline for the list of source names.

    Parameters
    ----------
    source_names_filename : str
        Filename for source names.
    path : str
        Path to save the output.
    """
    src_names = pd.read_csv(source_names_filename, header=None)
    run_for_sources(src_names[0], path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Mutually exclusive group for filename and source name
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--filename', help='Name of the file to read names of sources.')
    group.add_argument('--source_name', help='Name of the source.')
    group.add_argument('--ra', type=float, help='Right ascension of the source (in degrees).', required=False)

    # Path argument is not required, and declination should be provided in case usage --ra argument.
    parser.add_argument('--dec', type=float, help='Declination of the source (in degrees).', required=False)
    parser.add_argument('--path', type=str, help='declination of the source (in degrees)', required=False, default='./')

    args = parser.parse_args()
    source_name = args.source_name
    ra = args.ra
    dec = args.dec
    filename = args.filename
    path = args.path

    if source_name:
        run_for_sources([source_name], path=path)
    elif filename:
        read_source_names_from_file(filename, path=path)
    elif ra:
        if not dec:
            print('Please also provide the declination of the source.')
        else:
            run_with_coordinates(SkyCoord(f"{ra} {dec}", unit="deg"), path=path)
    else:
        print("You passed false arguments.")
