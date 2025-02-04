import os
import pandas as pd
from pyopenms import *

class IntensityCollector:
    def __init__(self, mzml_dir, nist_results_file, excluded_ions=None):
        self.mzml_dir = mzml_dir
        self.nist_results_file = nist_results_file
        self.excluded_ions = excluded_ions if excluded_ions is not None else []
        self.intensity_data = []  # To store the final intensity data

    def process_directory(self):
        """Process all mzML files in the specified directory."""
        mzml_files = [f for f in os.listdir(self.mzml_dir) if f.endswith(".mzML")]

        # Load NIST results once and structure it for easier lookup
        nist_matches = self.load_nist_results()

        # Initialize data structure to store output
        excel_data = {"File": []}  # The first column will be the file names

        for match in nist_matches:
            compound = match["name"]
            excel_data[compound] = {
                "Match Score": match["match_score"],
                "Ion": None,  # To be updated per file
                "Intensities": []
            }

        # Process each mzML file
        for mzml_file in mzml_files:
            mzml_path = os.path.join(self.mzml_dir, mzml_file)
            print(f"Processing: {mzml_file}")

            # Add file name to the list
            excel_data["File"].append(os.path.basename(mzml_file))

            # Extract intensities
            self.process_mzml_file(mzml_path, excel_data, nist_matches)

        # Convert to DataFrame and save as an Excel file
        self.save_to_excel(excel_data)

    def process_mzml_file(self, mzml_file, excel_data, nist_matches):
        """Process a single mzML file and collect intensities based on NIST search results."""
        exp = MSExperiment()
        MzMLFile().load(mzml_file, exp)

        for match in nist_matches:
            compound = match["name"]
            retention_times = match["retention_times"]
            top_ions = match["top_ions"]

            best_ion = self.get_best_ion(top_ions)
            if best_ion:
                mz, intensity = best_ion

                # Store ion only once (for first file processed)
                if excel_data[compound]["Ion"] is None:
                    excel_data[compound]["Ion"] = mz

                # Append intensity
                excel_data[compound]["Intensities"].append(intensity)
            else:
                # If no valid ion found, store NaN
                excel_data[compound]["Intensities"].append(np.nan)

    def save_to_excel(self, excel_data):
        """Save the collected intensities to an Excel file with the desired structure."""
        # Convert structured dictionary into a DataFrame
        files = excel_data["File"]

        # Prepare DataFrame rows
        df_rows = []
        df_rows.append(
            ["Match Score"] + [excel_data[compound]["Match Score"] for compound in excel_data if compound != "File"])
        df_rows.append(["Ion"] + [excel_data[compound]["Ion"] for compound in excel_data if compound != "File"])

        for i, file in enumerate(files):
            row = [file] + [excel_data[compound]["Intensities"][i] for compound in excel_data if compound != "File"]
            df_rows.append(row)

        # Convert to DataFrame
        df = pd.DataFrame(df_rows)

        # Save to Excel
        output_file = os.path.join(self.mzml_dir, "intensity_results.xlsx")
        df.to_excel(output_file, index=False, header=False)

        print(f"Intensity data saved to: {output_file}")

    def get_best_ion(self, top_ions):
        """Return the most abundant ion from the top ions list, excluding those on the exclusion list."""
        for ion in top_ions:
            mz, intensity = ion
            if mz not in self.excluded_ions:
                return ion
        return None  # Return None if no valid ion is found

    def load_nist_results(self):
        """Load and parse the NIST search results from the text file."""
        matches = []
        with open(self.nist_results_file, "r") as file:
            lines = file.readlines()

        current_match = None
        for line in lines:
            if line.startswith("Match:"):
                if current_match:
                    matches.append(current_match)
                current_match = {"name": line.split(":")[1].strip(), "retention_times": [], "top_ions": []}
            elif line.startswith("Retention Times:"):
                rt_str = line.split(":")[1].strip()
                retention_times = [float(rt) for rt in rt_str.strip("[]").split(",")]
                current_match["retention_times"] = retention_times
            elif line.startswith("Top Ions:"):
                ions_str = line.split(":")[1].strip()
                top_ions = [tuple(map(float, ion.split(","))) for ion in ions_str.strip("[]").split("),")]
                current_match["top_ions"] = top_ions

        if current_match:
            matches.append(current_match)

        return matches
