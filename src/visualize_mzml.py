import matplotlib.pyplot as plt
import numpy as np
from mzml_processor import MzMLProcessor
from peakutils import baseline
from scipy.signal import savgol_filter

mzml_data = r'C:\Jack\Projects\Metabolomics Project\11_11_24 TMS Std Curve\SC3 TMS.mzml'

int_matrix = MzMLProcessor.create_intensity_matrix(mzml_data)

scan_times = [float(meta["scan_start_time"]) for meta in int_matrix.spectra_metadata]
total_ion_currents = [float(meta["total_ion_current"]) for meta in int_matrix.spectra_metadata]
tic_array = np.array(total_ion_currents)

chunk_size = len(tic_array) // 2
num_chunks = len(tic_array) // chunk_size
if len(tic_array) % chunk_size != 0:
    num_chunks += 1

# Array to store estimated baseline values
baseline_values = np.zeros_like(tic_array)

# Minimum chunk size
min_chunk_size = chunk_size // 2

# Degrees for the function used for baseline correction
deg = 3

# Iterate over the chunks
for i in range(num_chunks):
    start_idx = i * chunk_size
    end_idx = (i + 1) * chunk_size

    # Ensure the last chunk doesn't go out of bounds
    if end_idx > len(tic_array):
        end_idx = len(tic_array)

    # Intensities for this chunk
    segment_intensities = tic_array[start_idx:end_idx]

    # If this is the last chunk and it's smaller than the threshold, combine it with the previous chunk
    if i == num_chunks - 1 and len(segment_intensities) < min_chunk_size:
        # Combine with the second-to-last chunk
        prev_chunk_end_idx = start_idx
        prev_chunk_start_idx = (i - 1) * chunk_size
        if prev_chunk_end_idx < len(tic_array):  # If not the very first chunk
            segment_intensities = np.concatenate([
                tic_array[prev_chunk_start_idx:prev_chunk_end_idx],
                segment_intensities
            ])

            # Calculate the baseline for the combined chunk
            bl = baseline(segment_intensities, deg=deg)

            # Store the baseline values
            baseline_values[prev_chunk_start_idx:end_idx] = bl
            continue  # Skip processing the last chunk as it's already combined

    # Apply baseline correction to the chunk
    bl = baseline(segment_intensities, deg=deg)

    # Store the baseline values
    baseline_values[start_idx:end_idx] = bl

smoothed_baseline = savgol_filter(baseline_values,window_length = 25, polyorder = 3)

plt.figure(figsize=(16,8))

plt.plot(scan_times, smoothed_baseline, color = 'green', linestyle='--', linewidth=2, label="Baseline")
plt.plot(scan_times, total_ion_currents, color='blue', linewidth=1.5)

plt.xlabel("Time")
plt.ylabel("Intensity")
plt.ylim(0,100000)
plt.title("Total Ion Current")
plt.grid(True)
plt.show()