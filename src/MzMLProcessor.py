"""Class that will accept a .d GCMS directory asn an input and convert it to .mzML for use in analysis"""
import subprocess
import os

class MzMLProcessor:
    def __init__(self, file_path):
        """Initialize with a .d or .ms file path"""
        self.file_path = file_path
        self.mzml_path = None  # Will be set after conversion

    """Convert .d or .ms file to .mzML using ProteoWizard's msconvert"""
    """outputs the .mzML file to the .d directory that the .mzML is made from"""
    def convert_to_mzml(self):

        # Set the output directory to be the same as the input file's directory
        output_dir = os.path.dirname(self.file_path)

        self.mzml_path = os.path.join(output_dir, os.path.basename(self.file_path) + ".mzML")

        command = [
            "msconvert", self.file_path, "--mzML", "--outdir", output_dir
        ]

        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"Conversion successful: {self.mzml_path}")
        else:
            print(f"Error converting file: {result.stderr}")
