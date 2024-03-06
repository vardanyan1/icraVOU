import os
import re
import sys
import errno
import logging
import argparse
import subprocess

import pandas as pd


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class Config:
    RADIUS = "4.0"
    EXCLUDED_CATALOGS = ["DEBL", "FMonLC"]


class FileManager:
    def check_file_exists(self, path):
        return os.path.isfile(path)

    def move_file(self, src, dest):
        os.rename(src, dest)

    def extract_value_from_wisestat(self, filepath):
        with open(filepath, 'r') as f:
            content = f.read()

        if "Warning" in content:
            return "N/A"

        match = re.search(r'Log\(Nu_peak\)\s*:\s*([\d\.\+\-]+) \+/- ([\d\.]+)', content)
        return "{} +/- {}".format(match.group(1), match.group(2)) if match else "N/A"


class CSVHandler:
    def concatenate_csv_files_with_pandas(self, files, output_file):
        df_list = [pd.read_csv(file) for file in files if not pd.read_csv(file).empty]
        combined_df = pd.concat(df_list, ignore_index=True)
        self.drop_invalid_rows(combined_df)
        combined_df.to_csv(output_file, index=False)

    def drop_invalid_rows(self, df):
        # Define columns to check for numeric values
        columns_to_check = ['freq. ', 'flux ', 'err_flux ']  # Ensure columns are correctly named as in your dataframe

        for column in columns_to_check:
            # Attempt to convert each value to numeric, invalid entries become NaN
            original_column = df[column].copy()
            df[column] = pd.to_numeric(df[column], errors='coerce')

            # Identify rows where conversion failed
            failed_conversion = df[df[column].isna() & original_column.notna()]
            if not failed_conversion.empty:
                for index, row in failed_conversion.iterrows():
                    message = "Invalid value in catalog '{}' for column '{}': {}".format(row['catalog'], column,
                                                                                         row[column])
                    ErrorHandler.log_info(message)

            # Drop rows with NaN values in the specified columns
            df.dropna(subset=[column], inplace=True)

        return df

    def handle_csv_data(self, filepath, ra, dec, radius):
        # Load the CSV file
        df = pd.read_csv(filepath, encoding='utf-8')

        # Sort the DataFrame by frequency
        df = df.sort_values(by='freq. ')

        # Exclude specific catalogs if needed
        for catalog in Config.EXCLUDED_CATALOGS:
            df = df[~df['catalog'].str.contains(catalog, case=False, na=False)].copy()

        # Format specified columns as scientific notation
        columns_to_format = ['freq. ', 'flux ', 'err_flux ']
        for column in columns_to_format:
            # Ensure the column exists in the DataFrame to avoid KeyErrors
            if column in df.columns:
                df[column] = df[column].apply(lambda x: "{:.4e}".format(x))

        # Define the output path for the processed CSV file
        # Updated to include RA, DEC, and RADIUS in the filename
        output_filename = "/work_dir/{}_{}_{}.csv".format(ra, dec, radius)

        # Write the processed DataFrame to a CSV file
        df.to_csv(output_filename, index=False)

        logging.info("Processed CSV data has been written to {}".format(output_filename))


class CommandExecutor:
    def execute(self, cmd):
        process = subprocess.Popen(cmd)
        process.wait()

    def change_working_directory(self, path):
        os.chdir(path)


class ErrorHandler:
    @staticmethod
    def log_info(message):
        logging.info(message)

    @staticmethod
    def log_error(message):
        logging.error(message)

    @staticmethod
    def log_warning(message):
        logging.warning(message)


class VouBlazarsHandler:

    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec
        self.cmd_executor = CommandExecutor()
        self.file_manager = FileManager()
        self.csv_handler = CSVHandler()
        self.error_handler = ErrorHandler()

        # Define the attributes here with default values
        self.w_peak = None
        self.csv_filepath = None

    def process(self):
        self.execute_vou_blazars_generation()
        self.handle_files()
        self.execute_scripts()
        self.handle_csv_tasks()

    def execute_vou_blazars_generation(self):
        self.error_handler.log_info("Executing command for VOU-Blazars generation...")
        cmd = ["./bin/vou-blazars-hybrid.sh", "--ra", self.ra, "--dec", self.dec, "--mode", "-s", "--allcats", "y"]
        self.cmd_executor.execute(cmd)
        self.error_handler.log_info("Command executed successfully.")

    def handle_files(self):

        self.csv_filepath = os.path.join("/VOU_Blazars/tmp", "Sed.csv")
        src_txt_filepath = os.path.join("/VOU_Blazars/tmp", "Sed.txt")

        if not self.file_manager.check_file_exists(src_txt_filepath):
            self.error_handler.log_error("File not found: {}".format(src_txt_filepath))
            raise IOError(errno.ENOENT, os.strerror(errno.ENOENT), src_txt_filepath)

        self.error_handler.log_info("handle files Done...")

    def execute_data_queries(self):
        """
        Execute the conesearch.py script to query both Fermi and UVOT data.
        """
        self.error_handler.log_info("Executing Fermi and UVOT data queries...")

        # Change directory to conesearch_files
        self.cmd_executor.change_working_directory("/VOU_Blazars/conesearch_files")

        # Updated command to include both Fermi and UVOT data paths
        cmd = ["python3", "conesearch.py",
               "--ra", str(self.ra), "--dec", str(self.dec),
               "--radius", Config.RADIUS]

        self.cmd_executor.execute(cmd)
        output_filename = os.path.join("/VOU_Blazars/conesearch_files", "{}_{}_{}.csv".format(
            self.ra, self.dec, Config.RADIUS))

        # Move the file back to the main directory
        self.file_manager.move_file(output_filename, os.path.join("../", output_filename))

        # Change directory back to the main directory
        self.cmd_executor.change_working_directory("..")

        self.error_handler.log_info("Data queries executed successfully.")
        return output_filename

    def execute_scripts(self):

        csv_filepath = os.path.join("/VOU_Blazars/tmp", "Sed.csv")
        if not self.file_manager.check_file_exists(csv_filepath):
            self.error_handler.log_error("File not found: {}".format(csv_filepath))
            raise IOError(errno.ENOENT, os.strerror(errno.ENOENT), csv_filepath)

        if float(self.dec) > -31:
            self.error_handler.log_info("Executing Zwicky data loader script...")
            self.cmd_executor.execute(["python3", "./zwicky_data_loader.py", "--ra", self.ra, "--dec", self.dec])
            self.error_handler.log_info("Zwicky data loader script executed.")

    def handle_csv_tasks(self):
        ztf_csv_filepath = "./ZTF_data.csv"
        combined_csv_filepath = "./combined_data.csv"
        # Execute the data queries and get the output filename
        combined_data_filepath = self.execute_data_queries()

        # Initialize the list of files to concatenate with Sed.csv which always exists
        files_to_concatenate = [self.csv_filepath]

        # Check and add ZTF_data.csv if it exists
        if self.file_manager.check_file_exists(ztf_csv_filepath):
            files_to_concatenate.append(ztf_csv_filepath)
        else:
            self.error_handler.log_warning("ZTF_data.csv not available for processing.")

        # Check and add combined_data.csv if it exists
        if self.file_manager.check_file_exists(combined_data_filepath):
            files_to_concatenate.append(combined_data_filepath)
        else:
            self.error_handler.log_warning("combined_data.csv not available for processing.")

        # If there are additional files to concatenate with Sed.csv
        if len(files_to_concatenate) > 1:
            self.error_handler.log_info("Handling combined data...")
            self.csv_handler.concatenate_csv_files_with_pandas(files_to_concatenate, combined_csv_filepath)
            self.csv_handler.handle_csv_data(combined_csv_filepath, self.ra, self.dec, Config.RADIUS)
            self.error_handler.log_info("Finished handling combined data.")
        else:
            # Only Sed.csv is available; handle it directly
            self.error_handler.log_info("Handling Sed.csv...")
            self.csv_handler.handle_csv_data(self.csv_filepath, self.ra, self.dec, Config.RADIUS)
            self.error_handler.log_info("Finished handling Sed.csv.")


try:  # wrap the whole script in a try-except block

    if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("--ra", type=str, required=True, help="RA value")
        parser.add_argument("--dec", type=str, required=True, help="Dec value")
        args = parser.parse_args()

        vou_handler = VouBlazarsHandler(args.ra, args.dec)
        vou_handler.process()

except IOError as e:
    # If the error is a file not found error
    if e.errno == errno.ENOENT:
        sys.stderr.write("FileNotFoundError occurred: {0}\n".format(e))
        sys.exit(2)  # return 2 to indicate a file not found error
    else:
        # Some other kind of IOError occurred
        sys.stderr.write("An IOError occurred: {0}\n".format(e))
        sys.exit(1)

except Exception as e:
    # Catching all other exceptions
    sys.stderr.write("An error occurred: {0}\n".format(e))
    sys.exit(1)
