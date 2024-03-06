import errno
import argparse
import sys
import os
import glob
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u


def filter_rows(df, target_ra, target_dec, max_distance):
    # Create SkyCoord objects for the entire DataFrame
    df_coords = SkyCoord(ra=df['Ra'].values * u.degree, dec=df['Dec'].values * u.degree)
    target_coord = SkyCoord(ra=target_ra * u.degree, dec=target_dec * u.degree)

    # Calculate the angular separation between the target and all rows
    separation = target_coord.separation(df_coords)

    # Filter based on max distance
    return df[separation < max_distance * u.arcsec]


def query_data(ra, dec, radius, path, catalog_name):
    data_df = pd.read_csv(path)
    filtered_df = filter_rows(data_df, ra, dec, radius)

    if not filtered_df.empty:
        # Use .copy() to ensure we're working with a DataFrame copy and avoid SettingWithCopyWarning
        filtered_df = filtered_df.copy()

        # Use .loc to safely assign new columns
        filtered_df.loc[:, 'catalog'] = catalog_name
        filtered_df.loc[:, 'reference'] = catalog_name

        # Use .rename method without inplace to avoid potential issues and return the modified DataFrame
        filtered_df = filtered_df.rename(columns={
            'Frequency': 'freq. ',
            'Flux': 'flux ',
            'dFlux': 'err_flux ',
            'MJD_start': 'MJD_start ',
            'MJD_end': 'MJD_end ',
            'Flag': 'flag'
        })

    return filtered_df


def process_all_csv(ra, dec, radius, directory):
    csv_files = glob.glob(os.path.join(directory, '*.csv'))
    dfs = []

    for file_path in csv_files:
        catalog_name = os.path.basename(file_path).split('.')[0]  # Assuming file name is catalog name
        dfs.append(query_data(ra, dec, radius, file_path, catalog_name))

    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df


def main(ra, dec, radius):
    directory = os.path.dirname(os.path.realpath(__file__))
    combined_df = process_all_csv(ra, dec, radius, directory)
    output_filename = f"{ra}_{dec}_{radius}.csv"

    if not combined_df.empty:
        output_columns = ['freq. ', 'flux ', 'err_flux ', 'MJD_start ', 'MJD_end ', 'flag', 'catalog', 'reference']
        existing_columns = [col for col in output_columns if col in combined_df.columns]
        combined_df.to_csv(output_filename, index=False, columns=existing_columns)

    return output_filename


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query astronomical data from CSV files.")
    parser.add_argument("--ra", type=float, required=True, help="RA value in degrees")
    parser.add_argument("--dec", type=float, required=True, help="Dec value in degrees")
    parser.add_argument("--radius", type=float, default=4.0, help="Radius value in arcsec")

    args = parser.parse_args()

    try:
        output_file = main(args.ra, args.dec, args.radius)
        print(f"Data combined and output to {output_file}")
    except IOError as e:
        if e.errno == errno.ENOENT:
            sys.stderr.write(f"FileNotFoundError occurred: {e}\n")
            sys.exit(2)
        else:
            sys.stderr.write(f"An IOError occurred: {e}\n")
            sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"An error occurred: {e}\n")
        sys.exit(1)
