import subprocess
from pathlib import Path
import numpy as np
import xml.etree.ElementTree as ET
import base64
import zlib
from intensity_matrix import IntensityMatrix

class MzMLProcessor:

    """Initialize with a .d file path"""
    def __init__(self, file_path):

        self.file_path = Path(file_path)  # Convert to Path object for easy manipulation
        self.mzml_path = file_path  # Will be set after conversion, set to file_path temporarily for testing (set back to None in final version)

    """Getters and setters"""
    #region
    @property
    def file_path(self):
        return self._file_path

    @file_path.setter
    def file_path(self,value):

        path = Path(value)

        if not path.exists():
            raise ValueError(f"File path does not exist: {value}")
        if not path.is_dir():
            raise ValueError(f"Expedcted a .d directory, but got a file: {value}")
        else:
            self._file_path = value
    #endregion

    """Convert .d or .ms file to .mzML using ProteoWizard's msconvert"""
    def convert_to_mzml(self):

        parent_dir = self.file_path.parent.parent  # Moves up one level from input file directory
        output_dir = parent_dir / "mzML_files"  # Creates a new "mzML_files" directory

        output_dir.mkdir(exist_ok=True)

        self.mzml_path = output_dir / (self.file_path.stem + ".mzML")  # Keep the same filename, change extension

        command = [
            "msconvert", str(self.file_path), "--mzML", "--outdir", str(output_dir)
        ]

        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"Conversion successful: {self.mzml_path}")
        else:
            print(f"Error converting file: {result.stderr}")

    """Decodes base64, decompresses zlib, and converts to a NumPy array with an associated m/z and time lists."""
    @staticmethod
    def decode_binary_data(encoded_data, dtype):

        decoded = base64.b64decode(encoded_data)
        decompressed = zlib.decompress(decoded)

        return np.frombuffer(decompressed, dtype=dtype)

    """Binns masses -0.3 to +0.7 of integer values"""
    @staticmethod
    def bin_masses(unique_mzs, intensity_matrix):
        
        # create a list to hold binned_mz values
        binned_mzs = []

        # determine int bin values such that n is added to binned_mzs if n-0.3 < mz <= n+0.7
        for mz in unique_mzs:
            if (mz + 0.3) > (int(mz) + 1) and int(mz+1) not in binned_mzs:
                binned_mzs.append(int(mz)+1)
            if (mz - 0.7) <= int(mz) and int(mz) not in binned_mzs:
                binned_mzs.append(int(mz))
        
        # dictionary that will hold the binned_mz : indices of unique_mzs to bin
        bin_tracker = {}

        # intialize bin_tracker with empty lists for each binned_mz value
        for mz in binned_mzs:
            bin_tracker[mz] = []
        
        # iterate over unique_mzs
        for idx, unique_mz in enumerate(unique_mzs):
            # find correct bin for unique_mz value and add the index of that value to bin_tracker
            for binned_mz in binned_mzs:
                if (binned_mz - 0.3) < unique_mz <= (binned_mz + 0.7):
                    bin_tracker[binned_mz].append(idx)

        _ , num_cols = np.shape(intensity_matrix)

        # create a new intensity matrix to hold binned values using the same number of columns but now with rows corrsponding to binned_mzs instead of unique_mzs
        binned_matrix = np.zeros((len(binned_mzs), int(num_cols)))

        # iterate over bin_tracker and sum rows of intesity_matrix as directed
        for row_idx, (binned_mz, indices_to_bin) in enumerate(bin_tracker.items()):

            # sum intensiteis of intensity_matrix for each time point
            summed_intensity = np.sum(intensity_matrix[indices_to_bin], axis=0)

            # add summed intensities to new matrix
            binned_matrix[row_idx] = summed_intensity

        return(binned_mzs, binned_matrix)

    """Extracts spectra metadata (scan start time, scan number, TIC value) and builds a matrix where each spectrum is
    represented by a column and each unique m/z is represented by a row.  If a spectrum does not include a value for an 
    m/z then it will be stored as a 0."""
    @staticmethod
    def create_intensity_matrix(mzml_path):

        tree = ET.parse(mzml_path)
        root = tree.getroot()

        namespaces = {
            '': 'http://psi.hupo.org/ms/mzml'
        }

        spectra_metadata = []
        intensity_list = []
        unique_mzs = set()

        # Iterate over each <spectrum> element
        for spectrum in root.findall('.//spectrum', namespaces):
            scan_id = spectrum.get('id')
            if scan_id:
                scan_id = scan_id.split('=')[-1]  # Split by '=' and take the second part (which is the number)

            scan_start_time = None
            scan_list = spectrum.find('scanList', namespaces)
            if scan_list is not None:
                scan = scan_list.find('scan', namespaces)
                if scan is not None:
                    scan_start_time = scan.find('.//cvParam[@name="scan start time"]', namespaces)
                    if scan_start_time is not None:
                        scan_start_time = scan_start_time.get('value')

            binary_data_elements = spectrum.findall('./binaryDataArrayList/binaryDataArray/binary', namespaces)

            if len(binary_data_elements) < 2:
                print(f"Skipping spectrum {scan_id} due to missing binary data")
                continue

            # Save metadata for the spectrum in the list
            metadata = {
                'scan_id': scan_id,
                'scan_start_time': scan_start_time,
            }
            spectra_metadata.append(metadata)

            # Grabs the binary encoded m/z and intensity arrays
            mz_array_encoded = binary_data_elements[0].text
            intensity_array_encoded = binary_data_elements[1].text

            mz_array = MzMLProcessor.decode_binary_data(mz_array_encoded, dtype=np.float64)
            intensity_array = MzMLProcessor.decode_binary_data(intensity_array_encoded, dtype=np.float32)

            # Create a dictionary for each spectrum (m/z -> intensity)
            spectrum_intensity_dict = dict(zip(mz_array, intensity_array))

            # Add the spectrum dictionary to the intensity_list
            intensity_list.append(spectrum_intensity_dict)

            # Add the m/z values to the set of unique m/z values
            unique_mzs.update(mz_array)

        # Convert the set of unique m/z values to a sorted list
        unique_mzs = sorted(unique_mzs)

        # Create a NumPy array for the intensity matrix
        intensity_matrix = np.zeros((len(unique_mzs), len(intensity_list)))

        # Create a map from m/z values to row indices in the intensity matrix
        mz_index_map = {mz: idx for idx, mz in enumerate(unique_mzs)}

        # Iterate over the intensity_list to fill the intensity matrix
        for col_idx, spectrum_intensity_dict in enumerate(intensity_list):
            for mz, intensity in spectrum_intensity_dict.items():
                row_idx = mz_index_map.get(mz)
                if row_idx is not None:
                    intensity_matrix[row_idx, col_idx] = intensity

        # bin the intensity matrix and unique_mzs
        binned_mzs, binned_matrix = MzMLProcessor.bin_masses(unique_mzs, intensity_matrix)

        # add TIC row to end of matrix
        sum_row = np.sum(binned_matrix, axis=0)
        final_matrix = np.vstack((binned_matrix,sum_row))

        # add 9999 value to end of binned_mzs to represent the TIC
        binned_mzs.append(9999)

        # create intensity matrix object
        output_matrix = IntensityMatrix(final_matrix,binned_mzs,spectra_metadata)

        # calculate the abundance thereshold
        output_matrix.calculate_threshold()

        # replace any values at or below abundance_threshold with corrosponding abundance threshold value
        output_matrix.apply_threshold()

        # calculate the noise factor
        output_matrix.calculate_noise_factor()

        return output_matrix
        