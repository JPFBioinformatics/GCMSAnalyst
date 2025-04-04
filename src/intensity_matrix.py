import numpy as np
from peakutils import baseline
from scipy.signal import savgol_filter, find_peaks

# Class for storage and cleaning of intensity matrix extracted by mzml_processor
class IntensityMatrix:

    def __init__(self, intensity_matrix, unique_mzs, spectra_metadata = None):
        self.intensity_matrix = intensity_matrix
        self.unique_mzs = unique_mzs
        self.spectra_metadata = spectra_metadata
        self.noise_factor = None
        self.abundance_threshold = None
        self.peak_list = None

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
            raise ValueError(f"unique m/z length {len(value)} does not match intensity array row count {self.intensity_matrix.shape[0]}")
        if not isinstance(value, list):
            raise ValueError('unique m/z is not a list')
        else:
            self._unique_mzs = value

    @property
    def spectra_metadata(self):
        return self._spectra_metadata

    @spectra_metadata.setter
    def spectra_metadata(self,value):
        if value is not None:
            if not len(value) == self.intensity_matrix.shape[1]:
                raise ValueError('Spectra metadata length does not match intensity array column count')
            if not isinstance(value, list):
                raise ValueError('Spectra metadata is not a list')
        self._spectra_metadata = value
    #endregion

    # region Abundance Threshold

    # calculates At value to replace 0 values with
    def calculate_threshold(self):

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

        intensity_matrix = self.intensity_matrix

        # the minimum measured abundance in the intensity matrix
        min_value = np.min(intensity_matrix[intensity_matrix>0])

        # creates a 2d array to store the threshold transitions, min value in the input array/sqrt(fraction of transitions in segment in row)
        threshold_values = np.empty((len(self.unique_mzs),10))

        # split the array into 10 approximately equal time segments
        segments = np.array_split(intensity_matrix, 10, axis=1)

        # counter for the start index of each segment
        start_idx = 0

        # list to hold the start index of each segment
        segment_starts = []

        # for each segment count the number of times a 0 value is followed by a nonzero value and store in transitions array
        for seg_idx,segment in enumerate(segments):
            for row_idx, row in enumerate(segment):

                # create array with 1 in any position where a 0 to nonzero transition occurs
                transitions = ((row[:-1] == 0) & (row[1:] > 0))
                # number of 0 to nonzero transitions in row
                num_transitions = np.sum(transitions)
                # length of the segment
                segment_length = segment.shape[1]
                # fraction of all scans in segment that are involved in m/z transtion
                threshold_values[row_idx, seg_idx] = num_transitions / segment_length

            # adds the start index for this segment to segment_starts list
            segment_starts.append(start_idx)
            # increments start_idx so that it now holds the first index value of the next segment
            start_idx += segment_length

        # takes the square root of all transition fraction values
        threshold_values **= 0.5

        # multiplies these square rooted values by the min value in matrix
        threshold_values *= min_value

        # dictionary that holds the start index of each segment (list) and the 2D numpy array (10 col, len(unique_mzs) rows) with threshold values stored in each cell
        threshold_dict = {
            'start_idxs' : segment_starts,
            'values' : threshold_values
        }
 
        self.abundance_threshold = threshold_dict

    # takes any value in the array that is below At for that segment for that m/z value and 
    def apply_threshold(self):

        matrix = self.intensity_matrix

        # iterate over each row, segment
        for row_idx,row in enumerate(matrix):
            for seg_idx in range(10):
                
                # start index for this segment
                start = self.abundance_threshold['start_idxs'][seg_idx]

                # end index for this segment
                if seg_idx <9:
                    end = self.abundance_threshold['start_idxs'][seg_idx+1]
                else:
                    end = row.shape[0]

                # get threshold value for this row this segment
                threshold = self.abundance_threshold['values'][row_idx,seg_idx]

                # replace values in this range for this segment with threshold
                row[start:end] = np.where(row[start:end]<threshold,threshold,row[start:end])

        self.intensity_matrix = matrix 

    # endregion

    # region Noise Factor Calculation

    # calculates the noise factor (Nf) for the entire intensity_matrix
    def calculate_noise_factor(self):

        matrix = self.intensity_matrix

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
                crossings = self.count_crossings(row,avg)

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
                current_nf = self.calculate_row_nf(row)
                # append result to median_devs list
                segment_nfs.append(current_nf)

            # adds the segment noise factors to the master list
            noise_factors.append(segment_nfs)

            self.noise_factor = np.median(noise_factors)

            return f"Noise Factor calculated: Nf = {self.noise_factor}"
    
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

        # calculat noise factor
        nf = np.median(deviations)/sqrt_of_mean
        
        # return the median of the deviations / sqrt of the mean (Nf for that row)
        return nf

    # endregion

    # region Finding Maxima

    # finds the peaks (maxima and bounds) for each row of a given intensity matrix and the tic, last row is TIC
    def identify_peaks(self, matrix, prom=None):

        # array to hold the lists of peak values, each entry of peaks corrosponds to a single m/z row in same order as unique_mzs
        peaks = []

        for row in matrix:
            row_peaks = self.find_maxima(row,prom=prom)
            peaks.append(row_peaks)

        self.peak_list = peaks

        return peaks

    # finds local maxima and bounds of peaks for a given 1D array
    def find_maxima(self, array, prom = None):

        # set default prominance
        if prom == None:
            prom = self.noise_factor*100

        # Excludes the first and last 12 points from the search to prevent bounding errors
        range = array[12:-12]

        # finds the local maxima of the given array, stores their index
        max_idxs, _ = find_peaks(range, prominence=prom)

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
            
            # width filter, skip peak if peak has less than 3 scans on either side of the peak
            if right_bound - max < 3 or max - left_bound < 3:
                continue

            # calculate baseline
            ten_baseline = self.tentative_baseline(left_bound,right_bound,array)
            baseline = self.least_squares_baseline(left_bound,right_bound,array,ten_baseline)

            # calcluate quadratic fit for peak
            fit = self.quadratic_fit(array,max)

            # grab the precise location and height of the peak
            precise_max_location = fit['x_values'][1]
            precise_max_height = fit['y_values'][1] - (baseline['slope']*(precise_max_location-baseline['left_bound']) + baseline['y_int'])
            precise_max_abundance = fit['y_values'][1]
            
            # finds the bin (0.1 of a scan) that the precise max is located within by truncating at 1 decimal point
            max_bin = int(precise_max_location*10) / 10

            # calculate convolution value for this peak
            conv = self.convolution_value(array,max)
            
            max_info = {
                'left_bound' : left_bound,
                'right_bound' : right_bound,
                'center' : max,
                'precise_max_location' : precise_max_location,
                'precise_max_height' : precise_max_height,
                'max_abundance' : precise_max_abundance,
                'bin' : max_bin,
                'conv_value' : conv
            }

            # accept peak if it passes threshold check
            if self.threshold_check(array,max,precise_max_height):
                maxima.append(max_info)

        # returns list of dictionary entries containing maxima information
        return maxima

    # finds the left or right deconvolution bound for a given maxima, step = 1 for right bound step = -1 for left bound
    def find_bound(self, array, center, step):

        nf = self.noise_factor
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
    
    # finds a quadratic fit for a set of 3 points in an array
    def quadratic_fit(self, array, center):

        # x values for fit, the center index and its two direct neighbors
        x_points = np.array([center-1, center, center+1])
        # y values for fit, from row corrosponding to x values
        y_points = array[x_points]

        # perform quadratic numpy polyfit, returning coefficients in a,b,c for ax^2 + bx + c form
        coeffs = np.polyfit(x_points, y_points, 2)
        a,b,c = coeffs

        # calcluate the precise maxima of the fit
        max_x = -b/(2*a)

        # array of left point(max_x-1), max, and right point (max_x +1)
        x_values = np.array([max_x-1, max_x, max_x+1]).astype(float)
        # array of y values corrosponding to x_values 
        y_values = a*x_values**2 + b*x_values + c

        fit_result = {
            'x_values' : x_values,
            'y_values' : y_values,
            'coeffs' : coeffs
        }

        return fit_result
    
    # checks if peak is above rejection threshold (4 noise units is base but we can adjust in the future)
    def threshold_check(self, row, peak_idx, height):

        threshold = 4 * self.noise_factor * row[peak_idx]**0.5 

        if height < threshold:
            return False
        else:
            return True

    # calculates convolution value for a single peak, used to see if peak is singlet or not
    def convolution_value(self,row,max):
        
        # holds the sum of all rates of sharpness calculated for the peak
        rate_sum = 0

        # loop over all the scans in the 3 scan window and calculate rate for each
        for i in range(1,4):
            term1 = (row[max+(i+1)] - row[max+i]) / row[max+i]
            term2 = (row[max-(i+1)] - row[max-i]) / row[max-i]

            rate_sum += term1+term2

        return rate_sum

    # endregion

    # region Baseline Calculation

    # calculates a tentative baseline for a percieved component
    def tentative_baseline(self,left_bound,right_bound,array):

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

        # shfit baseline down if any of its values are greater than the value of input array at same index
        for idx, element in enumerate(tentative_baseline):
            if element > component_array[idx]:
                diff = element - component_array[idx]
                tentative_baseline -= diff

        return tentative_baseline

    # calculates a least squars baseline based on a component array and its tentative baseline
    def least_squares_baseline(self,left_bound,right_bound,array,ten_baseline):
        
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
        baseline_array = m * np.arange(len(component_array)) + b

        # create dictioary to hold values
        bl = {
            'left_bound' : left_bound,
            'right_bound' : right_bound,
            'baseline_array' : baseline_array,
            'slope' : m,
            'y_int' : b
        }
        
        return bl

    # endregion

    # region Component Perception

    # finds the bin (1 scan = 10 bins) to place maxima in based on sharpness of peak
    def bin_peaks(self):
        
        dict_list = self.peak_list

        # iterates over all peaks in a list of lists of dictionaries dict_list to find the true bin to put that peak in
        for row_idx, peaks in enumerate(dict_list):
            for peak in peaks:

                # value that stores how much the bin the pak is currently in should be shifted by
                bin_offset = 0

                # maximum sharpness value calculated within the range
                max_sharpness = self.adjacent_bin_sharpness(row_idx,peak)

                # range in which to calculate sharpness values to find correct bin
                range  = int(0.5+50/max_sharpness)

                # max abundance and location of peak 
                abuncance = peak['max_abuncance']
                location = peak['bin']

                # checks through all bins in range to see if any have a larger sharpness than the origonal bin
                for i in range(1,range+1):
                    right_sharpness = self.adjacent_bin_sharpness(row_idx,abuncance,location+0.1*i)
                    if right_sharpness > max_sharpness:
                        max_sharpness = right_sharpness
                        bin_offset = i
                
                    left_sharpness = self.adjacent_bin_sharpness(row_idx,abuncance,location-0.1*i)
                    if left_sharpness>max_sharpness:
                        max_sharpness = left_sharpness
                        bin_offset = -i

                # update peak dictionary parameters with the true bin and associated sharpness
                peak['bin'] += 0.1*bin_offset
                peak['sharpness'] = max_sharpness
                # number of bins (each bin is 0.1 of a scan) to look to for model peak generation
                peak['model_range'] = int(0.5+50/max_sharpness)
    
    # calculates the sharpness for a bin and the bins directly on its left and right and returns the sum
    def adjacent_bin_sharpness(self,row_idx,max_bin,max_location):
        left_sharpness = self.avg_sharpness(row_idx,max_bin,max_location,adjacency=(-0.1))
        right_sharpness = self.avg_sharpness(row_idx,max_bin,max_location,adjacency=(0.1))
        center_sharpness = self.avg_sharpness(row_idx,max_bin,max_location)

        return sum(left_sharpness,right_sharpness,center_sharpness)
    
    # calculates the average sharpness for a single maxima
    def avg_sharpness(self, row_idx, max_abundance,max_bin, step = 3, adjacency = 0):
        
        # intensity at max 
        A_max = max_abundance

        # find the bin we are working from
        max_bin += adjacency

        # row we are working on right now
        row = self.intensity_matrix[row_idx]

        l_sharpnesses = []
        r_shaprnesses = []

        for i in range(1,step+1):
            
            # x value that we need to genearte an estimate of intensity for
            right_x = max_bin+i
            left_x = max_bin-i

            # (x,y) paris for linear fit
            # right side
            right1_x = int(right_x)
            right1_y = row[right1_x]
            right2_x = right1_x+1
            right2_y = row[right2_x]
            # left side
            left1_x = int(left_x)
            left1_y = row[left1_x]
            left2_x = left1_x-1
            left2_y = row[left2_x]
            
            # generate fit
            r_m,r_b = self.linear_fit(right1_x,right1_y,right2_x,right2_y)
            l_m,l_b = self.linear_fit(left1_x,left1_y,left2_x,left2_y)

            # apply fit
            r_An = r_m*right_x + r_b
            l_An = l_m*left_x + l_b

            # calculate sharpness
            right = self.calculate_sharpness(A_max,r_An,i)
            left = self.calculate_sharpness(A_max, l_An,i)

            # append shaprnesses to lists
            r_shaprnesses.append(right)
            l_sharpnesses.append(left)
        
        r_max = max(r_shaprnesses)
        l_max = max(l_sharpnesses)

        return ((r_max + l_max)/2)
    
    # creates a linear fit based on 2 points
    def linear_fit(self,x1,y1,x2,y2):

        slope = (y2-y1)/(x2-x1)
        y_int = y1-slope*x1

        return slope,y_int
    
    # calculates the shaprness for a a point with A_n intensity n scans from the peak max with intensity A_max
    def calculate_sharpness(self,A_max,A_n,n):

        sharpness = (A_max-A_n)/(n*self.noise_factor*(A_max**0.5))

        return sharpness
    
    
    
    # endregion

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
