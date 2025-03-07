import numpy as np
import matplotlib.pyplot as plt
from src.intensity_matrix import IntensityMatrix
from src.mzml_processor import MzMLProcessor

def generate_synthetic_chromatogram(length=5000, num_peaks=300, noise_level=10):
    x = np.linspace(50, length, length)
    chromatogram = np.zeros_like(x)

    # Generate random peaks
    for _ in range(num_peaks):
        peak_center = np.random.randint(100, length - 100)  # Avoid edges
        peak_width = np.random.randint(5, 20)
        peak_height = np.random.randint(150, 10000)

        peak = peak_height * np.exp(-((x - peak_center) ** 2) / (2 * peak_width ** 2))
        chromatogram += peak

    # Add baseline drift (slow, smooth change)
    baseline_drift = 20 * np.sin(0.0005 * np.pi * x)
    chromatogram += baseline_drift

    # Add random noise
    noise = np.random.normal(0, noise_level, size=length)
    chromatogram += noise

    return x, chromatogram

"""
# Generate test chromatogram
x, chromatogram = generate_synthetic_chromatogram()
first_correction = IntensityMatrix.correct_baseline(chromatogram, deg = 8, chunks = 1, low_chunk = 150)
difference = chromatogram - first_correction
difference = np.maximum(chromatogram - first_correction, 0)


# Plot
plt.figure(figsize=(10, 5))
plt.xlabel("Time")
plt.ylabel("Intensity")
plt.title("Synthetic Chromatogram for Baseline correction testing")

plt.plot(x, chromatogram, label="Original Chromatogram", color = 'blue', alpha = 0.5)
plt.plot(x,first_correction, label="First Correction", color = 'green', alpha = 0.8)
plt.plot(x, difference, label = "Corrected Chromatogram", color = 'red', alpha = 0.5)
#plt.plot(x,second_correction, label="Second Correction", color = 'orange', alpha = 0.8 )
#plt.plot(x,third_correction, label="Third Correction", color = 'red',)

plt.legend()

plt.show()
"""

mzml_data = r'C:\Jack\Projects\Metabolomics Project\11_11_24 TMS Std Curve\SC3 TMS.mzml'

int_matrix = MzMLProcessor.create_intensity_matrix(mzml_data)

time = np.array([entry['scan_start_time'] for entry in int_matrix.spectra_metadata], dtype = float)
tic = np.array([entry['total_ion_current'] for entry in int_matrix.spectra_metadata], dtype = float)

# Region m/z array values WE NEED TO BIN RAW DATA
print(f"m/z list: {int_matrix.unique_mzs}")

# region Test to see how the length of time each spectrum accounts for
# Average of 0.1052 sec/scan +/- 7.815e-6 (0.074%)
spectra_width = []
for i in range(1,len(time)-1):
    duration = float(time[i]) - float(time[i-1])
    spectra_width.append(duration)

avg_time = np.mean(spectra_width)
avg_deviation = np.std(spectra_width)

print(f"Average duration of each scan: {avg_time}")
print(f"Standard deviation of average deviation: {avg_deviation}")
print(f"Percent deviation {round(avg_deviation*100/avg_time,3)}%")
# endregion

baseline = IntensityMatrix.correct_baseline(tic, deg = 8, chunks = 1, low_chunk = 150)

corrected_tic = tic - baseline
corrected_tic = np.maximum(corrected_tic, 0)

baseline2 = IntensityMatrix.correct_baseline(corrected_tic, deg = 2, chunks = 20, low_chunk = 150)

double_corrected_tic = corrected_tic - baseline2
double_corrected_tic = np.maximum(double_corrected_tic, 0)

fig,ax = plt.subplots(nrows = 3, ncols = 1, figsize = (10,10), sharex = True, sharey = True)

ax[0].plot(time, tic, label = "Origonal TIC", color = 'blue', alpha = 0.5)
ax[0].plot(time, baseline, label = "Baseline", color = 'green', alpha = 0.8 )
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Intensity')
ax[0].legend()

ax[1].plot(time, corrected_tic, label = "Corrected TIC", color = 'red', alpha = 0.5)
ax[1].plot(time, baseline2, label = "Second baseline", color = 'purple', alpha = 0.8)
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Intensity')
ax[1].legend()

ax[2].plot(time, double_corrected_tic, label = "Final TIC", color = 'orange', alpha = 0.8)
ax[2].set_xlabel('Time')
ax[2].set_ylabel('Intensity')
ax[2].legend()

fig.suptitle(f"SC3 TIC", fontsize=16, fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.95])

plt.show()

