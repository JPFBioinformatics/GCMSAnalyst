import os
from pyopenms import *

class NISTSearch:
    def __init__(self, mzml_file, nist_database, score_threshold=0.9):
        self.mzml_file = mzml_file
        self.nist_database = nist_database
        self.score_threshold = score_threshold
        # Determine the directory of the input mzML file
        self.directory = os.path.dirname(mzml_file)
        self.output_file = os.path.join(self.directory, "nist_search_results.txt")
        self.match_data = []

    def search_and_filter(self):
        """ Search the mzML file against the NIST database and filter results by score threshold. """
        # Load the mzML file containing experimental data
        exp = MSExperiment()
        MzMLFile().load(self.mzml_file, exp)

        # Peak picking (if not already done in previous steps)
        peak_picker = SimplePeakPicker()
        peak_picker.pickExperiment(exp)

        # Set up database search (using NIST or any other database)
        params = DatabaseSearcher()
        params.setDatabaseFile(self.nist_database)  # Path to your NIST database
        params.setTolerance(0.01)  # m/z tolerance (e.g., 0.01 Da)
        params.setIonizationMode("positive")  # Set ionization mode (positive or negative)
        params.setPrecursorMassTolerance(0.01)  # Tolerance for precursor mass (optional)

        # Perform the database search
        search_results = params.search(exp)

        # Filter matches based on the match score (using the threshold)
        filtered_results = [result for result in search_results if result.getScore() > self.score_threshold]

        # Store match data
        self.collect_match_data(filtered_results, exp)

        # Print out the filtered results
        print("Filtered Matches:")
        for result in filtered_results:
            print(f"Match: {result.getName()}, Score: {result.getScore()}")

        # Optionally, store filtered results for further analysis in a text file
        self.save_results(filtered_results)

        # Return filtered results
        return filtered_results

    def collect_match_data(self, filtered_results, exp):
        """ Collect retention times and top 10 abundant ions for each match. """
        for result in filtered_results:
            match_data = {
                "name": result.getName(),
                "score": result.getScore(),
                "retention_times": [],
                "top_ions": []
            }

            # Collect retention times from the experimental data
            for spectrum in exp.getSpectra():
                if result.getName() in spectrum.getMetaData().get("peptide_ref"):
                    # Assuming retention time is stored as "RT" in metadata
                    rt = spectrum.getRT()
                    match_data["retention_times"].append(rt)

                    # Get the top 10 most abundant ions in the spectrum
                    mz_values = spectrum.get_peaks()[0]
                    intensities = spectrum.get_peaks()[1]
                    top_ions = sorted(zip(mz_values, intensities), key=lambda x: x[1], reverse=True)[:10]
                    match_data["top_ions"] = top_ions

            self.match_data.append(match_data)

    def save_results(self, filtered_results):
        """ Save the filtered results and additional data to a text file. """
        with open(self.output_file, "w") as file:
            for result in filtered_results:
                file.write(f"Match: {result.getName()}, Score: {result.getScore()}\n")

                # Store additional data (retention times and top 10 ions)
                for match in self.match_data:
                    if match["name"] == result.getName():
                        file.write(f"Retention Times: {match['retention_times']}\n")
                        file.write(f"Top Ions: {match['top_ions']}\n")

        print(f"NIST search results saved to: {self.output_file}")
