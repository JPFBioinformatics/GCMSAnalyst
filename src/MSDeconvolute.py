import os
from pyopenms import *

class MSDeconvolute:
    def __init__(self, mzml_file):
        # Get the directory of the input mzML file
        self.mzml_file = mzml_file
        self.directory = os.path.dirname(mzml_file)
        self.deconvoluted_file = os.path.join(self.directory, "deconvoluted_output.mzML")

    def deconvolute(self):
        """Deconvolute the given mzML file using OpenMS and save the deconvoluted result to a new file in the same directory."""
        # Load the mzML file
        exp = MSExperiment()
        MzMLFile().load(self.mzml_file, exp)

        # Peak picking: Pick peaks from raw data (if not done before)
        peak_picker = SimplePeakPicker()
        peak_picker.pickExperiment(exp)

        # Perform the deconvolution (centroiding)
        deconvoluter = CentroidedSpectra()
        deconvoluter.deconvolute(exp)

        # Save the deconvoluted file in the same directory as the input file
        MzMLFile().store(self.deconvoluted_file, exp)
        print(f"Deconvoluted mzML saved as: {self.deconvoluted_file}")
