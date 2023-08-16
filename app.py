import argparse
import errno
import os
import sys
import subprocess
import shutil


def handle_voublazars(ra, dec):
    cmd = ["./bin/vou-blazars-hybrid.sh", "--ra", ra, "--dec", dec, "--mode", "-s"]
    process = subprocess.Popen(cmd)
    process.wait()

    ra_main, ra_fraction = ra.split('.')
    dec_main, dec_fraction = dec.split('.')

    if ra_main.startswith('-'):
        ra_main = 'm' + ra_main[1:]  # Skip the "-" sign

    if dec_main.startswith('-'):
        dec_main = 'm' + dec_main[1:]  # Skip the "-" sign

    # folder_name = "{}_{}_{}_{}_".format(ra_main, ra_fraction, dec_main, dec_fraction)

    csv_filename = "Sed.csv"
    csv_filepath = os.path.join("/VOU_Blazars/tmp", csv_filename)

    txt_filename = "Sed.txt"
    txt_filepath = os.path.join("/VOU_Blazars/tmp", txt_filename)

    # Define the new CSV and TXT filenames
    new_csv_filename = "{}.{}_{}.{}.csv".format(ra_main, ra_fraction, dec_main, dec_fraction)
    new_txt_filename = "{}.{}_{}.{}.txt".format(ra_main, ra_fraction, dec_main, dec_fraction)

    # Define the new CSV and TXT file paths in /work_dir
    new_csv_filepath = os.path.join("/work_dir", new_csv_filename)
    new_txt_filepath = os.path.join("/work_dir", new_txt_filename)

    # Copy and rename the CSV and TXT files to /work_dir
    if os.path.isfile(csv_filepath):
        shutil.copy2(csv_filepath, new_csv_filepath)
        print("CSV File copied and renamed to " + new_csv_filepath)
    else:
        print("CSV File does not exist: " + csv_filepath)

    if os.path.isfile(txt_filepath):
        shutil.copy2(txt_filepath, new_txt_filepath)
        print("TXT File copied and renamed to " + new_txt_filepath)
    else:
        print("TXT File does not exist: " + txt_filepath)


try:  # wrap the whole script in a try-except block

    if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("--ra", type=str, required=True, help="RA value")
        parser.add_argument("--dec", type=str, required=True, help="Dec value")
        args = parser.parse_args()

        handle_voublazars(args.ra, args.dec)

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
