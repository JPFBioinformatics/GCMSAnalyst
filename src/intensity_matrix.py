import numpy as np
from peakutils import baseline
from scipy.signal import savgol_filter, find_peaks

# Class for storage and cleaning of intensity matrix extracted by mzml_processor
class IntensityMatrix:

    def __init__(self, intensity_matrix, unique_mzs, spectra_metadata):
        self.intensity_matrix = intensity_matrix
        self.unique_mzs = unique_mzs
        self.spectra_metadata = spectra_metadata
        self.noise_factor = None
        self.abundance_threshold = None

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

    # region OG Baseline, probably gonna delete
    
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

    # region Baseline Helpers
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

    # endregion

    # calculates At value to replace 0 values with
    def calculate_threshold(self, intensity_matrix):

        """
        Counts the number of zero to nonzero transtions for each m/z in 10 approximately equally sized time segments then takes the square root
        of these values and multiplies it by the minimum abundance measured in the entire intensity matrix, use this value to replace 0 values

        Parameters:
            intensity_matrix (np.ndarray): 2D numpy array where each row corrosponds to a m/z chromatogram intensity profile and each column
                                           corrosponds to a scan
        Returns:
            threshold_values (np.ndarray): 2D numpy array with 10 columns (1 per segment) and a row for each unique m/z in the input matrix
                                           each entry corrosponds to the calculated threshold value for that m/z in that segment
        """

        # the minimum measured abundance in the intensity matrix
        min_value = np.min(intensity_matrix[intensity_matrix>0])

        # creates a 2d array to store the threshold transitions, min value in the input array/sqrt(fraction of transitions in segment in row)
        threshold_values = np.empty((len(self.unique_mzs),10))

        # split the array into 10 approximately equal time segments
        segments = np.array_split(intensity_matrix, 10, axis=1)

        # for each segment count the number of ties a 0 value is followed by a nonzero value and store in transitions array
        for seg_idx,segment in enumerate(segments):
            for row_idx, row in enumerate(segment):
                transitions = (row[:-1] == 0 & row[1:] > 0)
                # number of 0 to nonzero transitions in row
                num_transitions = np.sum(transitions)
                # length of the segment
                segment_length = segment.shape[1]
                # fraction of all scans in segment that are involved in m/z transtion
                threshold_values[row_idx, seg_idx] = num_transitions / segment_length

        # takes the squar root of all transition fraction values
        threshold_values **= 0.5

        # multiplies these square rooted values by the min value in matrix, resulting in At matrix
        threshold_values *= min_value

        return threshold_values

    # region Noise Factor Calculation

    # calculates the noise factor (Nf) for the entire intensity_matrix
    def calculate_noise_factor(self, intensity_matrix):
        
        # adds a row to be the TIC, sum of all m/z values for each scan
        sum_row = np.sum(intensity_matrix, axis=0)
        matrix = np.vstack(intensity_matrix,sum_row)

        # determines how many 13 scan segments that we will have, if the last segment is not full it is excluded from calculations
        num_segments = matrix.shape[1] // 13

        # create an empty list of segments, each of which will be a numpy array with a row for each m/z chromatogram and a column for each 13 scan segment
        segments = []

        #loop over the number of segments creating each segment in segments as we go
        for i in range(num_segments):
            start = i*13
            end = (i+1)*13
            segment = matrix[:, start:end]
            
            # filters out any rows that contain 0 values
            removed_zeros = segment[~np.any(segment == 0, axis = 1)]
            
            # stores the segments that have sufficient number of "crossings"
            crossing_filtered = []

            # filters out any rows that "cross" the average intensity 6 or fewer times
            for row in removed_zeros:
                avg = np.mean(row)
                crossings = IntensityMatrix.count_crossings(row,avg)

                if crossings > 6:
                    crossing_filtered.append(row)
                
            segments.append(np.array(crossing_filtered))
        
        # list to hold all the noise factors for each row of each segment 
        noise_factors = []

        # iterate through each segment
        for segment in segments:
            # stores the median deviation for current segment
            segment_nfs = []

            # iterate through all rows of the segment
            for row in segment:
                # calculate noise factor for current row
                current_nf = IntensityMatrix.calculate_row_nf(row)
                # append result to median_devs list
                segment_nfs.append(current_nf)

            # adds the segment noise factors to the master list
            noise_factors.append(segment_nfs)

            return np.median(noise_factors)
    
    # counts the number of times the values of an array "cross" a given average value
    def count_crossings(self,row,avg):
        crossings = 0
        for i in range(len(row)-1):
            if (row[i] < avg and row[i+1] > avg) or (row[i] > avg and row[i+1] < avg):
                crossings += 1
        return crossings

    # calculates and returns the median deviation for a given 1D array
    def calculate_row_nf(self, row):

        # calculate the mean of the row
        mean = np.mean(row)
        sqrt_of_mean = mean ** 0.5

        # calculate deviation from the mean for all members of row
        deviations = np.abs(row-mean)

        # return the median of the deviations / sqrt of the mean (Nf for that row)
        return np.median(deviations)/sqrt_of_mean

    # endregion

    # region Finding Maxima

    # finds the peaks (maxima and bounds) for each row of a given intensity matrix and the tic
    def find_peaks(self, intensity_matrix):
        
        # adds a TIC row
        sum_row = np.sum(intensity_matrix, axis=0)
        matrix = np.vstack(intensity_matrix,sum_row)

        # array to hold the lists of peak values, one list for each m/z row 
        peaks = np.empty_like(matrix.shape[0])

        for idx,row in enumerate(matrix):

            row_peaks = self.find_maxima(row)
            peaks[idx] = row_peaks
        
        return peaks

    # finds local maxima and bounds of peaks for a given 1D array
    def find_maxima(self, array):

        # Excludes the first and last 12 points from the search to prevent bounding errors
        range = array[12:-12]

        # finds the local maxima of the given array, stores their index
        max_idxs, _ = find_peaks(range, height = self.noise_factor * 5)

        # Shifts indices found in the range for use in the original array
        max_idxs += 12
        
        # list to hold dictionary entries containing left_bound, right_bound and center for each maxima
        maxima = []

        # go through each maxima in list, find its deconvolution window and check if sinal is high enough to be included
        for max in max_idxs:

            # find the left bound of the deconvolution window
            left_bound = self.find_bound(array,max,-1)
            # find the right bound of the deconvolution window
            right_bound = self.find_bound(array,max,1)

            max_bounds = {
                'left_bound' : left_bound,
                'right_bound' : right_bound,
                'center' : max
            }

            maxima.append(max_bounds)

        # returns list of dictionary entries containing left bound, right bound, and center for each max (index values)
        return maxima

    # finds the left or right deconvolution bound for a given maxima, step = 1 for right bound step = -1 for left bound
    def find_bound(self, array, center, step):

        nf = self.calculate_noise_factor(array)
        counter = 1 * step
        min_value = array[center]
        max_value = array[center]

        # iterate up to 12 setps in given direction from center
        while counter <= 12 and counter >= -12:
            
            # if the value at this step is less than the current min, set the min to this value
            if array[center+counter] < min_value:
                min_value = array[center+counter]

            # if the value at this step is less than 5% of close window here
            if array[center+counter] < 0.01*max_value:
                return center+counter
            
            # if the value at this step is more than 5 nf greater than the minimum close the window at the previous step
            if array[center+counter] > 5*nf+min_value:
                return center+counter-step
            
            # increment counter
            counter += step

        # if no previous checks returned a value close window at 12 steps from the max
        return center+step*12
    
    #def reject_peak(self, array, center, left_bound, right_bound):

    # endregion

    # finds a quadratic fit for a set of 3 points, returns a dictionary 3 x_values, 3 y_values (based on quadratic fit) and the a, b, c coefficients for the fit
    def quadratic_fit(self, row_idx, center_idxs):

        # pulls the row we are working on from the intensity matrix
        row = self.intensity_matrix[row_idx]

        # array for holding the results
        results = []

        # loops over all indices in center_idxs and calcuates their quadratic fit
        for center in center_idxs:
            
            # x values for fit, the center index and its two direct neighbors
            x_points = np.array([center-1, center, center+1])
            # y values for fit, from row corrosponding to x values
            y_points = row[x_points]

            # perform quadratic numpy polyfit, returning coefficients in a,b,c for ax^2 + bx + c form
            coeffs = np.polyfit(x_points, y_points, 2)
            a,b,c = coeffs

            # calcluate the precise maxima of the fit
            max_x = -b/(2*a)
            # array of left point(max_x-1), max, and right point (max_x +1)
            x_values = np.array([max_x-1, max_x, max_x+1])
            # array of y values corrosponding to x_values 
            y_values = a*x_values**2 + b*x_values + c

            fit_result = {
                'x_values' : x_values,
                'y_values' : y_values,
                'coeffs' : coeffs
            }
            
            results.append(fit_result)

        return results

    # calculates the uniqueness value for each m/z in 10 approximately equal time chunks
    def calculate_uniqueness(self, intensity_matrix):
        
        # sores the uniqueness values for each m/z for each segment
        uniqueness_values = np.zeros((len(self.unique_mzs),10))
        
        # splits array into 10 approximately equal segments
        segments = np.array_split(intensity_matrix)

        for seg_idx,segment in enumerate(segments):
            for row_idx, row in enumerate(segment):

                # counts nonzero elements in this row
                counter = np.count_nonzero(row)
                # length of the row
                row_size = len(row)

                # calculates the fraction of entries in row that are nonzero
                nonzero_fraction = counter/row_size

                # appends nonzero_fraction to proper place in the uniqueness_values matrix
                uniqueness_values[row_idx,seg_idx] = nonzero_fraction
        
        return uniqueness_values

    # region Baseline Calculation

    # calculates a tentative baseline for a percieved component
    def tentative_baseline(left_bound,right_bound,array):

        # create componenet array 
        component_array = array[left_bound:right_bound+1]

        # creates an x-values array to use later for baseline computing, each x value is just an index value for input array
        x = np.arange(len(component_array))

        # get the index of the peak maximum
        max_idx = np.argmax(component_array)

        # get the index values of the minimum on the left and on the right of the max
        left_idx = np.argmin(component_array[:max_idx])
        right_idx = np.argmin(component_array[max_idx:]) + max_idx

        # get the intensity values associated with both these minimums
        left_val = component_array[left_idx]
        right_val = component_array[right_idx]

        # get linear baseline variables
        m = (right_val - left_val) / (right_idx - left_idx) 
        b = left_val - m * left_idx

        # generate tentative baseline array
        tentative_baseline = m * x + b

        return tentative_baseline

    # calculates a least squars baseline based on a component array and its tentative baseline
    def least_squares_baseline(left_bound,right_bound,array, ten_baseline):
        
        # create componenet array 
        component_array = array[left_bound:right_bound+1]

        # recalculates abundances as height above tentative baseline
        adjusted_abundances = component_array - ten_baseline

        # get the indicies for the smallest 50% of values in the array
        sorted_indices = np.argsort(adjusted_abundances)
        smallest_indices = sorted_indices[:len(component_array)//2]

        # get y values (abundances) for least squares baseline estimation
        smallest_abundances = component_array[smallest_indices]
        
        # linear least squares fitting 
        m, b = np.polyfit(smallest_indices, smallest_abundances, 1)

        # generate a bseline for the entire compoonent array
        baseline = m * np.arange(len(component_array)) + b
        
        return baseline

    # endregion

    # calculates the sharpness for a single set of points
    def calculate_sharpness(self, row, start_idx, end_idx):
        
        sharpness = (row[start_idx] - row[end_idx]) / (abs(start_idx - end_idx) * self.noise_factor * (row[start_idx]**0.5))

        return sharpness
