""" Takes .msp files from NIST database, extracted with Lib2NIST, and stores relevant data in JSON file for easy updat to SQL DB
    in Lib2NIST we binned masses with options -> add the following term to all m/z before rounding : 0.3 """

from pathlib import Path
from datetime import datetime
import json

class DBBuilder:

    def __init__(self,directory):
        """ Saves directory path as Path obj for easy utilization """
        self.directory = Path(directory)

    def parse_msp(self,file):
        """ Takes .msp file and converts each spectra to an entry in a json file """

        # create lists/dict to store records as we go
        records = []
        current_record = {}
        current_mzs = []
        current_intensities = []

        with open(file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()

                # if there is an empty line then the spectra data is done and move on to next spectra
                if not line:

                    # add the m/z list and intensity lists to current_record
                    if current_mzs and current_intensities:
                        current_record["mzs"] = current_mzs
                        current_record["intensities"] = current_intensities

                    # append current spectra info to records list
                    records.append(current_record)

                    # reset current_record peak info lists
                    current_record = {}
                    current_mzs = []
                    current_intensities = []
                    continue

                # grab the data from key:value pairs
                if ':' in line:

                    # identify key,value pairs
                    key,value = line.split(':',1)

                    # strip keys and values to get clean string
                    key = key.strip()
                    value = value.strip()

                    # store key value parirs
                    current_record[key] = value

                # if no ':' present then it is the m/z intensity pairs
                else:

                    # split the line at the space
                    mz,intensity = line.split(' ',1)

                    # store the m/z and intesnty values
                    current_mzs.append(int(mz))
                    current_intensities.append(int(intensity))

        # if there is no empty line after last record stor it anyway before we move on
        if current_record:
                if current_mzs and current_intensities:
                    current_record["mzs"] = current_mzs
                    current_record["intensities"] = current_intensities

        # return records list
        return records

    def parse_sdf(self,file):
        """ Takes .sdf file and converts each spectra to an entry in a json file """

        # create lists/dict to store records as we go
        records = []
        current_record = {}
        current_mzs = []
        current_intensities = []

        # open and decode file
        with open(file, 'r', encoding='utf-8') as f:

            # stores line number for debugging
            line_num = 1

            # reads current and next line into memory simultaneously
            line = next(f,None)
            next_line = next(f,None)

            while line is not None:
                line = line.strip()

                # if the line is "$$$$" this means this entry is done and to store data and move on
                if line == "$$$$":

                    # add the m/z list and intensity lists to current_record
                    if current_mzs and current_intensities:
                        current_record["mzs"] = current_mzs
                        current_record["intensities"] = current_intensities

                    # append current spectra info to records list
                    records.append(current_record)

                    # reset current_record peak info lists
                    current_record = {}
                    current_mzs = []
                    current_intensities = []

                    # advance loop
                    line = next_line
                    next_line = next(f,None)
                    line_num += 1
                    continue

                # grab the data from key:value pairs
                elif line.startswith('>'):

                    # identify key,value pairs
                    parts = line.split('  ',1)

                    # find and store a key
                    if len(parts) > 1:

                        # find key and remove <> surrounding it
                        key = parts[1].strip()
                        key = key[1:-1]

                    # for the mz/intensity pair data do a special parse to get values
                    if key == "MASS SPECTRAL PEAKS":
                        print("extracting ms info")
                        
                        # advance the loop
                        line = next_line
                        next_line = next(f,None)
                        line_num += 1

                        # loop over all m/z values:
                        while line is not None:
                            line = line.strip()

                            # check if there is some sort of error in file structure (unlikely as files are highly structured)
                            if line == ""  or line.startswith('>') or line == "$$$$":
                                print(f"SFD parse error: m/z, intensity paris in source file {file} corrupted at line {line_num}")
                                break
                            
                            # if structure looks fine continue the data gathering
                            # split line at the space
                            mz,intensity = line.split(' ',1)

                            # store m/z, intensity pair
                            current_mzs.append(int(mz))
                            current_intensities.append(int(intensity))

                            # advance loop 
                            line = next_line
                            next_line = next(f,None)
                            line_num += 1

                        # go to next line and continue outer loop
                        line = next_line
                        next_line = next(f,None)
                        line_num += 1
                        continue

                    # grab the values for the rest of the fields
                    elif next_line is not None:

                        # grab and store values
                        value = next_line.strip()
                        current_record[key] = value

                        # advance loop
                        line = next_line
                        next_line = next(f,None)
                        line_num += 1
                        continue

        # return records list
        return records

    def combine_to_json(self,output_json,type):
        """ Goes over all files in directory, removes data, and stores it all in one large json file
            can do msp or sdf, but sdf is default as it supplies both CasNo and INCHIKEY """

        # grab directory and start a list for all records to combine
        directory = self.directory
        all_records = []

        # for every file that ends in .msp first parse then add it to the all_records list
        for file in directory.glob(f"*.{type}"):

            # print message to ensure we can see if it fails
            print(f"Parsing {file}...")

            # save records for msp or sdf file
            if type == "msp":
                record = self.parse_msp(file)
            elif type == "sdf":
                record = self.parse_sdf(file)

            print(f"file {file} parsed")

            # add record to records
            all_records.extend(record)

        # change output to Path object
        output_path = Path(output_json)

        # open output and dump data to it
        with output_path.open('w',encoding='utf-8') as out:
            json.dump(all_records,out,indent=2)

        # print confirmation message
        print(f"Saved {len(all_records)} spectra to {output_path}")

if __name__ == "__main__":

    #get datetime when data extraction is started
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")

    # SDF files for my purposes
    type = "msp"

    # input directory
    dir = r"C:\Jack\code\Data\NIST2020_subset.MSP"

    # output path and name
    output = Path(r"C:\Jack\code\Data") / f"{timestamp}_{type}_NIST20.json"

    # build DB and ouptut json results
    builder = DBBuilder(dir)
    builder.combine_to_json(output,type)