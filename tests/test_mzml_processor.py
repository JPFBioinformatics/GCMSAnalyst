from pathlib import Path
from mzml_processor import MzMLProcessor
from intensity_matrix import IntensityMatrix

# region Creating an intensity matrix from mzML file

# get the root directory (GCMSAnalyst)
root_dir =  Path(__file__).resolve().parent.parent

# get the path to the example data file SC1_TMS.mzML
mzml_data = root_dir / 'example_data' / 'SC1_TMS.mzML'

# create an intensity matrix from the given mzML file
int_matrix = MzMLProcessor.create_intensity_matrix(mzml_data)
print(f"intensity matrix dimensions: {int_matrix.intensity_matrix.shape}")

