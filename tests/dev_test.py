import numpy as np
import matplotlib.pyplot as plt
from intensity_matrix import IntensityMatrix

# generates a synthetic chromatogram for testing
def generate_synthetic_chromatogram(length=3600, num_peaks=1, noise_level=500, peak_center = None, peak_width = None, peak_height = None):

    # generat x values and empty y values array
    x = np.arange(length)
    chromatogram = np.zeros_like(x).astype(np.float64)

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

    return x, chromatogram.reshape(1,-1)

# region Graphing synthetic chromatograms
"""
# generate synthetic chromatograms
x1, chromatogram_1 = generate_synthetic_chromatogram(peak_center=2000,peak_width=2,peak_height=500000)
x2, chromatogram_2 = generate_synthetic_chromatogram(peak_center=2006,peak_width=2,peak_height=400000)
combined_chromatogram = chromatogram_1 + chromatogram_2

# generate IntensityMatrix objects for each chromatogram
chrom1 = IntensityMatrix(chromatogram_1,[9999])
chrom2 = IntensityMatrix(chromatogram_2,[9999])
c_chrom = IntensityMatrix(combined_chromatogram,[9999])
print(f"chrom1 shape: {chrom1.intensity_matrix.shape}")
print(f"chrom2 shape: {chrom2.intensity_matrix.shape}")
print(f"c_chrom shape: {c_chrom.intensity_matrix.shape}")

# generate noise factors and abundance thresholds for chromatograms
chrom1.calculate_noise_factor()
chrom2.calculate_noise_factor()
c_chrom.calculate_noise_factor()
print(f"Chrom 1 Nf: {chrom1.noise_factor}")
print(f"Chrom 2 Nf: {chrom2.noise_factor}")
print(f"Convoluted chrom Nf: {c_chrom.noise_factor}")

# finds maxima and bounds for component array for chromatogram 1 and 2
c1_maxima = chrom1.identify_peaks(chrom1.intensity_matrix,prom=chrom1.noise_factor*1000)
c2_maxima = chrom2.identify_peaks(chrom2.intensity_matrix,prom=chrom2.noise_factor*1000)
conv_maxima = c_chrom.identify_peaks(c_chrom.intensity_matrix,prom=c_chrom.noise_factor*1000)
print("")
print("Maxima info:")
print(c1_maxima)
print(c2_maxima)
print(conv_maxima)
print("")

# calculate x values for compponent arrays used in graphing baselines
c1_x = np.arange(c1_maxima[0][0]['left_bound'], c1_maxima[0][0]['right_bound']+1)
c2_x = np.arange(c2_maxima[0][0]['left_bound'], c2_maxima[0][0]['right_bound']+1)
conv_x1 = np.arange(conv_maxima[0][0]['left_bound'],conv_maxima[0][0]['right_bound']+1)
conv_x2 = np.arange(conv_maxima[0][1]['left_bound'],conv_maxima[0][1]['right_bound']+1)

# calculate tentative baseline for chromatogram 1 and 2
tent_c1 = chrom1.tentative_baseline(c1_maxima[0][0]['left_bound'], c1_maxima[0][0]['right_bound'],chrom1.intensity_matrix[0])
tent_c2 = chrom2.tentative_baseline(c2_maxima[0][0]['left_bound'], c2_maxima[0][0]['right_bound'],chrom2.intensity_matrix[0])
tent_conv1 = c_chrom.tentative_baseline(conv_maxima[0][0]['left_bound'],conv_maxima[0][0]['right_bound'],c_chrom.intensity_matrix[0])
tent_conv2 = c_chrom.tentative_baseline(conv_maxima[0][1]['left_bound'],conv_maxima[0][1]['right_bound'],c_chrom.intensity_matrix[0])

# calculate least squares baseline for chromatogram 1 and 2
lsb_c1 = chrom1.least_squares_baseline(c1_maxima[0][0]['left_bound'], c1_maxima[0][0]['right_bound'],chrom1.intensity_matrix[0],tent_c1)
lsb_c2 = chrom2.least_squares_baseline(c2_maxima[0][0]['left_bound'], c2_maxima[0][0]['right_bound'],chrom2.intensity_matrix[0],tent_c2)
lsb_conv1 = c_chrom.least_squares_baseline(conv_maxima[0][0]['left_bound'],conv_maxima[0][0]['right_bound'],c_chrom.intensity_matrix[0],tent_conv1)
lsb_conv2 = c_chrom.least_squares_baseline(conv_maxima[0][1]['left_bound'],conv_maxima[0][1]['right_bound'],c_chrom.intensity_matrix[0],tent_conv2)

bl_center = (conv_maxima[0][0]['center'] - conv_maxima[0][0]['left_bound'])
peak1_height = c_chrom.intensity_matrix[0][conv_maxima[0][0]['center']] - lsb_conv1['baseline_array'][bl_center]

# calculate the height of each chromatogram and plot it as a green vertical line
chrom1_height = c1_maxima[0][0]['precise_max_height']
chrom1_location = c1_maxima[0][0]['precise_max_location']
chrom1_precise_max_baseline = lsb_c1['slope']*(chrom1_location-lsb_c1['left_bound'])+lsb_c1['y_int']
chrom1_precise_max_top = chrom1_precise_max_baseline + chrom1_height
print("Chromatogram 1 height information:")
print(f"Location: {chrom1_location}")
print(f"Height: {chrom1_height}")
print(f"Bottom of height line: {chrom1_precise_max_baseline}")
print(f"Top of height line: {chrom1_precise_max_top}")
print("")

chrom2_height = c2_maxima[0][0]['precise_max_height']
chrom2_location = c2_maxima[0][0]['precise_max_location']
chrom2_precise_max_baseline = lsb_c2['slope']*(chrom2_location-lsb_c2['left_bound'])+lsb_c2['y_int']
chrom2_precise_max_top = chrom2_precise_max_baseline + chrom2_height
print("Chromatogram 2 height information:")
print(f"Location: {chrom2_location}")
print(f"Height: {chrom2_height}")
print(f"Bottom of height line: {chrom2_precise_max_baseline}")
print(f"Top of height line: {chrom2_precise_max_top}")
print("")

conv1_height = conv_maxima[0][0]['precise_max_height']
conv1_location = conv_maxima[0][0]['precise_max_location']
conv1_precise_max_baseline = lsb_conv1['slope']*(conv1_location-lsb_conv1['left_bound'])+lsb_conv1['y_int']
conv1_precise_max_top = conv1_precise_max_baseline + conv1_height

conv2_height = conv_maxima[0][1]['precise_max_height']
conv2_location = conv_maxima[0][1]['precise_max_location']
conv2_precise_max_baseline = lsb_conv2['slope']*(conv2_location-lsb_conv2['left_bound'])+lsb_conv2['y_int']
conv2_precise_max_top = conv2_precise_max_baseline + conv2_height
print("Conv peaks height information:")
print(f"peak 1 height: {conv1_height}")
print(f"peak 2 height: {conv2_height}")

# Plot
fig,ax = plt.subplots(nrows = 3, ncols = 1, figsize = (10,10), sharex = True, sharey = True)

ax[0].plot(x1, c_chrom.intensity_matrix[0], label="Convoluted Chromatogram", color = 'blue', alpha = 0.5)
ax[0].plot(conv_x1, tent_conv1, label='Tentative baseline', color = 'grey', alpha = 0.5)
ax[0].plot(conv_x1, lsb_conv1['baseline_array'], label='Least Squares baseline', color = 'purple', alpha = 0.5)
ax[0].vlines(conv1_location,conv1_precise_max_baseline,conv1_precise_max_top,linestyles ='dashed',color = 'black')
ax[0].vlines(conv2_location,conv2_precise_max_baseline,conv2_precise_max_top,linestyles ='dashed',color = 'black')
ax[0].plot(conv_x2, tent_conv2, color = 'grey', alpha = 0.5)
ax[0].plot(conv_x2, lsb_conv2['baseline_array'], color = 'purple', alpha = 0.5)
ax[0].scatter(conv_maxima[0][0]['left_bound'], c_chrom.intensity_matrix[0][conv_maxima[0][0]['left_bound']], color = 'orange')
ax[0].scatter(conv_maxima[0][0]['right_bound'], c_chrom.intensity_matrix[0][conv_maxima[0][0]['right_bound']], color = 'orange')
ax[0].scatter(conv_maxima[0][0]['center'], c_chrom.intensity_matrix[0][conv_maxima[0][0]['center']], color = 'black')
ax[0].scatter(conv_maxima[0][1]['left_bound'], c_chrom.intensity_matrix[0][conv_maxima[0][1]['left_bound']], color = 'orange')
ax[0].scatter(conv_maxima[0][1]['right_bound'], c_chrom.intensity_matrix[0][conv_maxima[0][1]['right_bound']], color = 'orange')
ax[0].scatter(conv_maxima[0][1]['center'], c_chrom.intensity_matrix[0][conv_maxima[0][1]['center']], color = 'black')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Intensity')
ax[0].legend()


ax[1].plot(x1, chrom1.intensity_matrix[0], label="chromatogram 1", color = 'red', alpha = 0.5)
ax[1].plot(c1_x, tent_c1, label='Tentative baseline', color = 'grey', alpha = 0.5)
ax[1].plot(c1_x, lsb_c1['baseline_array'], label='Least Squares baseline', color = 'purple', alpha = 0.5)
ax[1].vlines(chrom1_location,chrom1_precise_max_baseline,chrom1_precise_max_top,linestyles ='dashed',color = 'black')
ax[1].scatter(c1_maxima[0][0]['left_bound'], chrom1.intensity_matrix[0][c1_maxima[0][0]['left_bound']], color = 'orange')
ax[1].scatter(c1_maxima[0][0]['right_bound'], chrom1.intensity_matrix[0][c1_maxima[0][0]['right_bound']], color = 'orange')
ax[1].scatter(c1_maxima[0][0]['center'], chrom1.intensity_matrix[0][c1_maxima[0][0]['center']], color = 'black')
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Intensity')
ax[1].legend()

ax[2].plot(x2, chrom2.intensity_matrix[0], label="chromatogram 2", color = 'green', alpha = 0.5)
ax[2].plot(c2_x, tent_c2, label='Tentative baseline', color = 'grey', alpha = 0.5)
ax[2].plot(c2_x, lsb_c2['baseline_array'], label='Least Squares baseline', color = 'purple', alpha = 0.5)
ax[2].vlines(chrom2_location,chrom2_precise_max_baseline,chrom2_precise_max_top,linestyles ='dashed',color = 'black')
ax[2].scatter(c2_maxima[0][0]['left_bound'], chrom2.intensity_matrix[0][c2_maxima[0][0]['left_bound']], color = 'orange')
ax[2].scatter(c2_maxima[0][0]['right_bound'], chrom2.intensity_matrix[0][c2_maxima[0][0]['right_bound']], color = 'orange')
ax[2].scatter(c2_maxima[0][0]['center'], chrom2.intensity_matrix[0][c2_maxima[0][0]['center']], color = 'black')
ax[2].set_xlabel('Time')
ax[2].set_ylabel('Intensity')
ax[2].legend()

plt.show()
"""
# endregion

