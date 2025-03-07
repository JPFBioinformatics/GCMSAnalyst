"""Class for storage and cleaning of intensity matrix extracted by mzml_processor"""

import numpy as np
from peakutils import baseline
from scipy.signal import savgol_filter


class IntensityMatrix:

    def __init__(self, intensity_matrix, unique_mzs, spectra_metadata):
        self.intensity_matrix = intensity_matrix
        self.unique_mzs = unique_mzs
        self.spectra_metadata = spectra_metadata

    """Getter and setter methods"""
    #region Getter/Setters
    @property
    def intensity_matrix(self):
        return self._intensity_matrix

    @intensity_matrix.setter
    def intensity_matrix(self,value):
        if isinstance(value,np.ndarray):
            self._intensity_matrix = value
        else:
            raise ValueError("intensity matrix is not a numpy array")

    @property
    def unique_mzs(self):
        return self._unique_mzs

    @unique_mzs.setter
    def unique_mzs(self,value):
        if not len(value) == self.intensity_matrix.shape[0]:
            raise ValueError('unique m/z length does not match intensity array row count')
        if not isinstance(value, list):
            raise ValueError('unique m/z is not a list')
        else:
            self._unique_mzs = value

    @property
    def spectra_metadata(self):
        return self._spectra_metadata

    @spectra_metadata.setter
    def spectra_metadata(self,value):
        if not len(value) == self.intensity_matrix.shape[1]:
            raise ValueError('Spectra metadata length does not match intensity array column count')
        if not isinstance(value, list):
            raise ValueError('Spectra metadata is not a list')
        else:
            self._spectra_metadata = value
    #endregion

    """Baseline correction"""
    # Uses correct_baseline method to iteratively correct the baseline noise for a given array, not the whole matrix
    @staticmethod
    def iterative_baseline_correction(intensity_array):

        # Compute first, single chunk baseline estimation
        smoothed_baseline = IntensityMatrix.correct_baseline(intensity_array, deg = 8, chunks = 1, low_chunk = 100)

        # Generates a corrected intensity array based on computed baseline
        corrected_array = intensity_array - smoothed_baseline
        # Sets any negative values in corrected array to 0
        corrected_array = np.maximum(corrected_array, 0)

        # Runs a second baseline correction with a finer chunk size and a lower degree polynomial fit
        second_baseline = IntensityMatrix.correct_baseline(corrected_array, deg = 2, chunks = 20)

        # Generates a final corrected array based on the resultant corrected array from first correction
        final_array = corrected_array - second_baseline
        # Sets any negative values in final array to 0
        final_array = np.maximum(final_array, 0)

        return final_array

    # region Helper Functions
    @staticmethod
    def correct_baseline(intensity_array, deg = 8, chunks = 1, low_chunk = 100):

        # Stores the baseline correction values array returned by peakutils baseline() function
        baseline_values = np.zeros_like(intensity_array)

        # Set the size of the chunks you will be processing
        chunk_size = len(intensity_array) // chunks
        # Ensures the chunk size is within reasonable bounds
        if chunk_size < low_chunk:
            chunk_size = low_chunk
        if chunk_size > 1500:
            chunk_size = 1500

        # Sets min_chunk_size, if a chunk is smaller than this it will be combined with previous chunk
        min_chunk_size = chunk_size // 2

        # Sets the number of chunks, ensuring an extra is added if needed
        num_chunks = len(intensity_array) // chunk_size
        if len(intensity_array) % chunk_size != 0:
            num_chunks += 1

        # Iterate over the chunks
        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = (i + 1) * chunk_size

            # Ensure the last chunk doesn't go out of bounds
            if end_idx > len(intensity_array):
                    end_idx = len(intensity_array)

            # Intensities for this chunk
            segment_intensities = intensity_array[start_idx:end_idx]

            # If this is the last chunk, and it's smaller than the threshold combine it with the previous chunk
            if i == num_chunks - 1 and len(segment_intensities) < min_chunk_size:
                # Combine with the second-to-last chunk
                prev_chunk_end_idx = start_idx
                prev_chunk_start_idx = (i - 1) * chunk_size
                if prev_chunk_end_idx < len(intensity_array):  # If not the very first chunk
                    segment_intensities = np.concatenate([
                        intensity_array[prev_chunk_start_idx:prev_chunk_end_idx],
                        segment_intensities
                    ])

                    # Calculate the baseline for the combined chunk
                    bl = baseline(segment_intensities, deg = deg)


                    # Store the baseline values
                    baseline_values[prev_chunk_start_idx:end_idx] = bl
                    #Removes negative values from basleine correction array
                    baseline_values = np.maximum(baseline_values,0)

                    continue  # Skip processing the last chunk as it's already combined

            # Apply baseline correction to the chunk
            chunk_baseline = baseline(segment_intensities, deg = deg)

            # Store the baseline values
            baseline_values[start_idx:end_idx] = chunk_baseline
            #removes negative values from baseline correction array
            baseline_values = np.maximum(baseline_values,0)

        # Uses Savitzky-Golay smoothing to make sure chunk edges merge smoothly
        sg_window = IntensityMatrix.sg_window_length(chunk_size)

        smoothed_baseline = savgol_filter(baseline_values, window_length = sg_window, polyorder = 2)

        return smoothed_baseline

    # calculates a dynamic window length for smoothing algorithm
    @staticmethod
    def sg_window_length(chunk_size):
        wl = chunk_size // 10
        if wl % 2 == 0:
             wl += 1
        return wl

    @staticmethod
    def polynomial_correct_negative_values(array, window_size = 5, degree = 1):

        # Creates a list of all indices of array where the value is negative
        negative_indices = np.where(array < 0)[0].astype(int)

        # Returns original array if no negative indices
        if len(negative_indices) == 0:
            return array

        clusters = []
        current_cluster = [negative_indices[0]]

        for i in range(1, len(negative_indices)):  # Iterate through negative indices starting from the second element
            if negative_indices[i] == negative_indices[i - 1] + 1:
                # If the current index is exactly 1 greater than the previous, it's part of the same cluster
                current_cluster.append(negative_indices[i])
            else:
                # If the current index is NOT consecutive, store the current cluster and start a new one
                clusters.append(current_cluster)
                current_cluster = [negative_indices[i]]  # Start a new cluster with the current index

        #adds the last cluster of negative values to the list
        clusters.append(current_cluster)

        # Iterates through each cluster of negative indices, taking window_size points from each side of the cluster,
        # then only uses positive values in that segment to create an n degree polynomial fit (linear, 1, is usually
        # plenty) for that segment and then replaces negative values in that segment with predictions from the fit
        for cluster in clusters:
            start_idx = max(int(cluster[0]) - window_size, 0)
            end_idx = max(int(cluster[-1] + 1) + window_size, len(array))

            segment = array[start_idx : end_idx]
            x = np.arange(start_idx, end_idx)
            y = segment

            valid_x = x[y >= 0]
            valid_y = y[y >= 0]

            if len(valid_x) > degree:
                fit = np.polyfit(valid_x, valid_y, deg = degree)
                fit_line = np.polyval(fit,x)

                array[start_idx : end_idx] = np.where(
                    array[start_idx : end_idx] < 0, fit_line, array[start_idx : end_idx]
                )

        return array
    #endregion

    """Chromatogram Smoothing"""
