import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from intensity_matrix import IntensityMatrix
from mzml_processor import MzMLProcessor
from pathlib import Path

# region Creating an intensity matrix from mzML file

# get the root directory (GCMSAnalyst)
root_dir =  Path(__file__).resolve().parent.parent

# get the path to the example data file SC1_TMS.mzML
mzml_data = root_dir / 'example_data' / 'SC1_TMS.mzML'

# create an intensity matrix from the given mzML file
int_matrix = MzMLProcessor.create_intensity_matrix(mzml_data)

# endregion

# region Extract data from the IntensityMatrix object

# Scan numbers
scan_no = np.array([entry['scan_id'] for entry in int_matrix.spectra_metadata], dtype = int)

# Time
time = np.array([entry['scan_start_time'] for entry in int_matrix.spectra_metadata], dtype = float)
print(f"Number of time points sampled: {len(time)}")

# Total ion current
tic = np.array([entry['total_ion_current'] for entry in int_matrix.spectra_metadata], dtype = float)

# Intensity matrix
intensity_matrix = int_matrix.intensity_matrix

# extract the list of m/z values from the IntensityMatrix object
mz_list = int_matrix.unique_mzs
print(f"Number of unique m/z values: {len(mz_list)}")

# endregion

# region Determine the time accounted for in each scan, average: 0.1052 sec/scan +/- 7.815e-6 (0.074%)

spectra_width = []
for i in range(1,len(time)-1):
    duration = float(time[i]) - float(time[i-1])
    spectra_width.append(duration)

avg_time = np.mean(spectra_width)
avg_deviation = np.std(spectra_width)

print(f"Average duration of each scan: {avg_time}")
print(f"Standard deviation of average deviation: {avg_deviation}")
print(f"Percent deviation {round(avg_deviation*100/avg_time,3)}%")

# endregion

# region Graphing the TIC

# create plot
plt.plot(scan_no,tic)

# set layout
plt.title("Total Ion Current")
plt.xlabel("Scan")
plt.ylabel("Intensity")

plt.show()

# endregion
