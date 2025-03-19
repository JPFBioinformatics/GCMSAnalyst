import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# calculates the noise factor for a given chromatogram
def calculate_noise_factor(array):
    
    # split array into 13 scan segments excluding the last segment if it has < 13 elements for simplicity
    num_segments = len(array) // 13

    # create a list to hold the segment sub-arrays
    segments = []

    # loop num_segments times to create subarrays
    for i in range(num_segments):
        start = i*13
        end = (i+1)*13
        segment = array[start:end]

        # calculate the average for the segment
        avg = np.mean(segment)

        # appends segment only if it "crosses" the average > 6 times
        if count_crossings(segment,avg) > 6:
            segments.append(segment) 

    # list to hold calculated noise factors for calculated segments
    noise_factors = []

    # iterates over each segment and calculates the noise factor
    for segment in segments:

        # calculates noise factor
        nf = calculate_row_nf(segment)

        # appends noise factor to the master list
        noise_factors.append(nf)

    return np.median(noise_factors)

# calculates and returns the median deviation for a given 1D array
def calculate_row_nf(row):

    # calculate the mean of the row
    mean = np.mean(row)
    sqrt_of_mean = mean ** 0.5

    # calculate deviation from the mean for all members of row
    deviations = np.abs(row-mean)

    # return the median of the deviations / sqrt of the mean (Nf for that row)
    return np.median(deviations)/sqrt_of_mean

# counts the number of times an array "crosses" a given average value
def count_crossings(row,avg):
    crossings = 0
    for i in range(len(row)-1):
        if (row[i] < avg and row[i+1] > avg) or (row[i] > avg and row[i+1] < avg):
            crossings += 1
    return crossings

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

# finds local maxima for a given 1D array
def find_maxima(array):

    # Excludes the first and last 12 points from the search to prevent bounding errors
    range = array[12:-12]

    # finds the local maxima of the given array, stores their index
    max_idxs, _ = find_peaks(range, height = 500)

    # Shifts indices found in the range for use in the original array
    max_idxs += 12
        
    # list to hold dictionary entries containing left_bound, right_bound and center for each maxima
    maxima = []

    # go through each maxima in list, find its deconvolution window and check if sinal is high enough to be included
    for max in max_idxs:

        # find the left bound of the deconvolution window
        left_bound = find_bound(array,max,-1)

        # find the right bound of the deconvolution window
        right_bound = find_bound(array,max,1)

        max_bounds = {
            'left_bound' : left_bound,
            'right_bound' : right_bound,
            'center' : max
        }

        maxima.append(max_bounds)

    # returns list of dictionary entries containing left bound, right bound, and center for each max (index values)
    return maxima

# finds the left or right deconvolution bound for a given maxima, step = 1 for right bound step = -1 for left bound
def find_bound(array, center, step):

    nf = calculate_noise_factor(array)
    counter = 1 * step
    min_value = array[center]
    max_value = array[center]

    # iterate up to 12 setps in given direction from center
    while counter <= 12 and counter >= -12:
            
        # if the value at this step is less than the current min, set the min to this value
        if array[center+counter] < min_value:
            min_value = array[center+counter]

        # if the value at this step is less than 5% of close window here
        if array[center+counter] < 0.05*max_value:
            return center+counter
            
        # if the value at this step is more than 5 nf greater than the minimum close the window at the previous step
        if array[center+counter] > 5*nf+min_value:
            return center+counter-step
            
        # increment counter
        counter += step

    # if no previous checks returned a value close window at 12 steps from the max
    return center+step*12

# generates a synthetic chromatogram for testing
def generate_synthetic_chromatogram(length=3600, num_peaks=1, noise_level=10, peak_center = None, peak_width = None, peak_height = None):

    # generat x values and empty y values array
    x = np.linspace(50, length, length)
    chromatogram = np.zeros_like(x)

    # Generate random peaks
    for _ in range(num_peaks):
        if peak_center is None:
            peak_center = np.random.randint(100, length - 100)  # Avoid edges
        if peak_width is None:
            peak_width = np.random.randint(5, 20)
        if peak_height is None:
            peak_height = np.random.randint(150, 10000)

        peak = peak_height * np.exp(-((x - peak_center) ** 2) / (2 * peak_width ** 2))
        chromatogram += peak

    # Add baseline drift (slow, smooth change)
    baseline_drift = 20 * np.sin(0.005 * np.pi * x)
    chromatogram += baseline_drift

    # Add random noise
    noise = np.random.normal(0, noise_level, size=length)
    chromatogram += noise

    # shift chromatogram up so that no values are negative
    min_value = np.min(chromatogram)
    if min_value < 0:
        chromatogram = chromatogram - min_value + 25
    if min_value > 0:
        chromatogram = chromatogram + 25

    return x, chromatogram

# generate synthetic chromatograms
x1, chromatogram_1 = generate_synthetic_chromatogram(peak_center=2000,peak_width=15,peak_height=7000)
x2, chromatogram_2 = generate_synthetic_chromatogram(peak_center=2015,peak_width=5,peak_height=1000)

# add chromatograms together to simulate convolution
combined_chromatogram = chromatogram_1 + chromatogram_2

# finds maxima and bounds for component array for chromatogram 1 and 2
c1_maxima = find_maxima(chromatogram_1)
c2_maxima = find_maxima(chromatogram_2)

print("chromatogram 1 bounds")
for key, value in c1_maxima[0].items():
    print(f"{key}: {value}")
print("chromatogram 2 bounds")
for key, value in c2_maxima[0].items():
    print(f"{key}: {value}")


# calculate x values for compponent arrays used in graphing baselines
c1_x = np.arange(c1_maxima[0]['left_bound'], c1_maxima[0]['right_bound']+1)
c2_x = np.arange(c2_maxima[0]['left_bound'], c2_maxima[0]['right_bound']+1)

# calculate tentative baseline for chromatogram 1 and 2
tent_c1 = tentative_baseline(c1_maxima[0]['left_bound'], c1_maxima[0]['right_bound'],chromatogram_1)
tent_c2 = tentative_baseline(c2_maxima[0]['left_bound'], c2_maxima[0]['right_bound'],chromatogram_2)

# calculate least squares baseline for chromatogram 1 and 2
lsb_c1 = least_squares_baseline(c1_maxima[0]['left_bound'], c1_maxima[0]['right_bound'],chromatogram_1,tent_c1)
lsb_c2 = least_squares_baseline(c2_maxima[0]['left_bound'], c2_maxima[0]['right_bound'],chromatogram_2,tent_c2)

# Plot
fig,ax = plt.subplots(nrows = 3, ncols = 1, figsize = (10,10), sharex = True, sharey = True)

ax[0].plot(x1, combined_chromatogram, label="Convoluted Chromatogram", color = 'blue', alpha = 0.5)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Intensity')
ax[0].legend()


ax[1].plot(x1, chromatogram_1, label="chromatogram 1", color = 'red', alpha = 0.5)
ax[1].plot(c1_x, tent_c1, label='Tentative baseline', color = 'grey', alpha = 0.5)
ax[1].plot(c1_x, lsb_c1, label='Least Squares baseline', color = 'purple', alpha = 0.5)
ax[1].scatter(c1_maxima[0]['left_bound'], chromatogram_1[c1_maxima[0]['left_bound']], color = 'black')
ax[1].scatter(c1_maxima[0]['right_bound'], chromatogram_1[c1_maxima[0]['right_bound']], color = 'black')
ax[1].scatter(c1_maxima[0]['center'], chromatogram_1[c1_maxima[0]['center']], color = 'black')
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Intensity')
ax[1].legend()

ax[2].plot(x2, chromatogram_2, label="chromatogram 2", color = 'green', alpha = 0.5)
ax[2].plot(c2_x, tent_c2, label='Tentative baseline', color = 'grey', alpha = 0.5)
ax[2].plot(c2_x, lsb_c2, label='Least Squares baseline', color = 'purple', alpha = 0.5)
ax[2].scatter(c2_maxima[0]['left_bound'], chromatogram_2[c2_maxima[0]['left_bound']], color = 'black')
ax[2].scatter(c2_maxima[0]['right_bound'], chromatogram_2[c2_maxima[0]['right_bound']], color = 'black')
ax[2].scatter(c2_maxima[0]['center'], chromatogram_2[c2_maxima[0]['center']], color = 'black')
ax[2].set_xlabel('Time')
ax[2].set_ylabel('Intensity')
ax[2].legend()

plt.show()
