import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from matplotlib import cm
from sympy.plotting.pygletplot.util import interpolate

# Constants and parameter
verbose = 2
register = np.arange(4095)

plot_spectra = True  # Bool whether spectra shall be created
test = True #Bool for testing, when true, only 1 file is analyzed

# Load flat file
print('Reading flat file ... ')
flatname = "C:/Users/Celina/Masterarbeit/Kartenmessung/Flattable/flat_241119_113313__absblack.tab_50.spec"
flatdata = np.genfromtxt(flatname)  # creates data for flattening with flatname file
print('Reading flatfile completed!')

#create flat file
mF = flatdata[:, 1]
norm = mF / np.median(mF[1704:1905])

# Manually create Mask
startpos = list(range(1, 84))  # prev 1

brightpos = list(range(1632, 1644 + 1))

triplepos = [103, 104, 105, 106, 107, 108, 109, 110, 111, 121, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140,
             157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195,
             231, 232, 244, 245, 256, 257, 258, 316, 317, 329, 330, 343, 401, 402, 414, 415, 428, 461, 486, 487, 488,
             499, 500, 501, 513, 514, 572, 573, 585, 586, 598, 599, 657, 658, 659, 669, 670, 671, 672, 674, 682, 683,
             684, 742, 743, 744, 755, 756, 768, 769, 801, 802, 803, 804, 826, 827, 828, 829, 839, 840, 841, 842, 852,
             853, 854, 855, 880, 881, 886, 887, 906, 907, 912, 913, 925, 926, 927, 938, 939, 940, 997, 998, 999, 1011,
             1012, 1024, 1025, 1083, 1084, 1096, 1097, 1109, 1110, 1168, 1169, 1181, 1182, 1194, 1195, 1235, 1248, 1266,
             1267, 1338, 1339, 1351, 1352, 1365, 1423, 1424, 1436, 1437, 1450, 1481, 1494, 1508, 1509, 1510, 1511, 1522,
             1523, 1535, 1536, 1593, 1594, 1595, 1596, 1606, 1607, 1608, 1609, 1619, 1620, 1621, 1622, 1645, 1646, 1647,
             1673, 1679, 1680, 1692, 1693, 1695, 1696, 1704, 1705, 1706, 1763, 1764, 1765, 1776, 1777, 1778, 1779, 1789,
             1790, 1791, 1792, 1796, 1809, 1823, 1849, 1850, 1862, 1863, 1875, 1876, 1934, 1935, 1947, 1948, 1961, 2019,
             2020, 2021, 2033, 2034, 2046, 2047, 2105, 2106, 2118, 2119, 2129, 2130, 2131, 2132, 2190, 2191, 2203, 2204,
             2216, 2217, 2287, 2288, 2289, 2290, 2300, 2301, 2302, 2303, 2315, 2361, 2374, 2387, 2445, 2446, 2458, 2459,
             2472, 2473, 2531, 2544, 2545, 2557, 2558, 2616, 2617, 2629, 2630, 2642, 2643, 2701, 2702, 2714, 2715, 2716,
             2727, 2728, 2786, 2787, 2799, 2800, 2812, 2813, 2871, 2872, 2884, 2885, 2898, 2957, 2969, 2970, 2983, 3042,
             3055, 3056, 3068, 3069, 3127, 3128, 3140, 3141, 3153, 3154, 3192, 3204, 3205, 3218, 3212, 3213, 3225, 3226,
             3227, 3238, 3239, 3297, 3298, 3299, 3309, 3310, 3311, 3312, 3323, 3324, 3467, 3468, 3469, 3480, 3481, 3482,
             3493, 3494, 3495, 3553, 3566, 3579, 3638, 3651, 3652, 3664, 3665, 3723, 3724, 3736, 3737, 3749, 3750, 3808,
             3809, 3821, 3822, 3834, 3835, 3894, 3907, 3920, 3978, 3979, 3991, 3992, 4005, 4015, 4028, 4041,
             4064]  # immer 3 Linien im Abstand 13

faintpos = [87, 94, 496, 539, 546, 615, 624, 625, 631, 709, 710, 736, 795, 871, 956, 965, 971, 972, 1050, 1051, 1057,
            1127, 1135, 1142, 1222, 1468, 1476, 1542, 1561, 1568, 1573, 1712, 1738, 1817, 1843, 1902, 1941, 1987, 1994,
            2072, 2073, 2101, 2151, 2156, 2157, 2184, 2185, 2212, 2223, 2249, 2275, 2276, 2328, 2419, 2420, 2425, 2498,
            2590, 2663, 2760, 2839, 2956, 3186, 3245, 3265, 3271, 3272, 3350, 3441, 3442, 3451, 3520, 3546, 3686, 3687,
            4077, 4090]  # oft 2er Paare

mask = startpos + brightpos + triplepos + faintpos
len_mask = len(mask)
#-----

def read_tables(path, table):

    with open(path + table, 'r') as file:

        ra_list = []
        dec_list = []
        spectra_list = []

        print('Start reading file ', table)
        for line in file:
            if line[0] == '#':
                continue
            if line[0] != 's':
                continue
            else:
                # print(line)
                spectrum_line = line[235:]
                spec = spectrum_line.split()
                if len(spec) < 4095:
                    continue
                else:
                    spectra_list.append(np.array(spec, dtype=float))
                    ra_list.append(float(line[120:127]))
                    dec_list.append(float(line[134:141]))

                del spec

        return ra_list, dec_list, spectra_list

def flat_spectra(spec_arr):
    norm_spec_arr = np.zeros(len(register))
    for i in register:
        if norm[i] > 0.:
            norm_spec_arr[i] = spec_arr[i] / norm[i]
        else:
            norm_spec_arr[i] = 0.
    return norm_spec_arr

def mask_spectra(spec_arr):
    spec_arr_masked = spec_arr.copy()
    for element in mask:
        spec_arr_masked[element + 1] = np.nan
    return spec_arr_masked


def interpolate_spectra(spec_arr):

    # Create a copy to avoid modifying the original array
    interpol_spec = spec_arr.copy()

    for idx in range(len(interpol_spec)):
        if np.isnan(interpol_spec[idx]):

            if idx == 0 and np.isnan(interpol_spec[idx]):
                interpol_spec[idx] = 0

            elif idx == len(interpol_spec) - 1 and np.isnan(interpol_spec[idx]):
                interpol_spec[idx] = interpol_spec[idx - 1]
            else:
                # Find the nearest previous non-NaN value
                left_idx = idx - 1
                # Find the nearest next non-NaN value
                right_idx = idx + 1
                while right_idx < len(interpol_spec) and np.isnan(interpol_spec[right_idx]):
                    right_idx += 1

                # Interpolate as the mean of neighboring values

                left_val = interpol_spec[left_idx]
                right_val = interpol_spec[right_idx]
                # Compute the mean of available neighbors
                interpol_spec[idx] = (left_val + right_val) / 2

    return interpol_spec


def norm_spectra(spec_arr):
    spec_norm = np.zeros(len(register))

    spec_median = np.median(spec_arr)

    counter = 0  # count nan's after normalization: if any element is nan, the median will be nan
    for i in register:
    # normalize spectra by their individual median
        spec_norm[i] = spec_arr[i] / spec_median

        if np.isnan(spec_norm[i]):
            counter += 1  # for debugging

    if counter != 0:  # for debugging
        print(n, ':', counter)  # for debugging

    return spec_norm

def sum_spectra(spec_arr, lower_lim, upper_lim, ra_list, dec_list, path, table):

    spectra_diff = spec_arr[lower_lim:upper_lim] - 1  # select area of data points and subtract 1
    spectra_sum = np.sum(spectra_diff)
    with open(path + table + '_sums.txt', 'a') as sum_table:
        sum_table.write("%s\t%s\t%s\n" % (ra_list[n], dec_list[n], spectra_sum))


# ------------ Start data analysis -----------------------
#path to project directory
project_path = "C:/Users/Celina/Masterarbeit/Kartenmessung/"

# Path to directory of map data tables
file_path = project_path + "Tabellen/"
#list of files in the directory
filelist = os.listdir(file_path)

#if test: reduce list to single file for debugging
if test:
    filelist = ['241108_112024___ALT39.5.tab']  # Used for analysis of a single file

#check length of filelist
len_filelist = len(filelist)
if verbose > 1:
    print(len(filelist))

#Loop for processing individual files
for tab in filelist:

    #read table
    RA_list, DEC_list, SPECTRA_list = read_tables(file_path, tab)
    print('Finished reading file: ', tab)

    #create plot figures
    fig1, ax1 = plt.subplots(1)  # figure for the raw spectra
    fig2, ax2 = plt.subplots(1)  # figure for the normalized spectra
    fig3, ax3 = plt.subplots(1)  # figure for the masked spectra
    fig4, ax4 = plt.subplots(1)  # figure for the interpolated spectra
    fig5, ax5 = plt.subplots(1)  # figure for the interpolated, individually normalized
    fig6, ax6 = plt.subplots(1)  # figure of a small range of the final spectra

    len_SPEC = len(SPECTRA_list)
    if verbose > 1:
        print('Total number of spectra in table:', len_SPEC)
        print('Control length of the spectrum (4096): ', len(SPECTRA_list[0]))

    k, j = 0, 0  # counter for every 1000 spectra to track progress
    n = 0  # counts number of spectra

    #loop for processing each spectrum in a table
    for spectrum in SPECTRA_list:

        k += 1  # Counter for number of spectra (for print only)
        if k == 1000:
            j += 1  # Counter for j * 1000 spectra, for print only
            print('Norm spectra', j * k)
            k = 0

        #flat spectrum
        flat_spec = flat_spectra(spectrum)

        #remove sharp peaks
        flat_spec_masked = mask_spectra(flat_spec)

        #interpolate spectra
        interpol = interpolate_spectra(flat_spec_masked)

        #normalize spectra by individual median
        flat_spec_masked_interpol_norm = norm_spectra(interpol)

        #sum the signal height in a specified range
        sum_path = project_path + "Sum_Tables/"
        sum_spectra(flat_spec_masked_interpol_norm, 1650,2000, RA_list, DEC_list, sum_path, tab)

        #plot several spectra after each processing step
        if plot_spectra:
            if 600 <= n <= 700:
                ax1.plot(spectrum, color='blue', alpha=0.01, lw=1)
                ax2.plot(flat_spec, color='blue', alpha=0.01, lw=1)
                ax3.plot(flat_spec_masked, color='blue', alpha=0.01, lw=1)
                ax4.plot(interpol, color='blue', alpha=0.01, lw=1)
                ax5.plot(flat_spec_masked_interpol_norm, color='blue', alpha=0.01, lw=1)
                ax6.plot(flat_spec_masked_interpol_norm, color='blue', alpha=0.01, lw=1)

            if 9600 <= n <= 9700:
                ax1.plot(spectrum, color='red', alpha=0.01, lw=1)
                ax2.plot(flat_spec, color='red', alpha=0.01, lw=1)
                ax3.plot(flat_spec_masked, color='red', alpha=0.01, lw=1)
                ax4.plot(interpol, color='red', alpha=0.01, lw=1)
                ax5.plot(flat_spec_masked_interpol_norm, color='red', alpha=0.01, lw=1)
                ax6.plot(flat_spec_masked_interpol_norm, color='red', alpha=0.01, lw=1)

            if 19600 <= n <= 19700:
                ax1.plot(spectrum, color='green', alpha=0.01, lw=1)
                ax2.plot(flat_spec, color='green', alpha=0.01, lw=1)
                ax3.plot(flat_spec_masked, color='green', alpha=0.01, lw=1)
                ax4.plot(interpol, color='green', alpha=0.01, lw=1)
                ax5.plot(flat_spec_masked_interpol_norm, color='green', alpha=0.01, lw=1)
                ax6.plot(flat_spec_masked_interpol_norm, color='green', alpha=0.01, lw=1)

        n += 1

    print('File processing complete:', tab)

    # Adjust plot aesthetics for all figures
    ax1.set_title('Raw Spectra')
    ax1.set_xlabel('Register Number')
    ax1.set_ylabel('Intensity')
    ax1.set_xlim((300, 3100))
    ax1.set_ylim((14000, 24000))

    ax2.set_title('Flattened Spectra')
    ax2.set_xlabel('Register Number')
    ax2.set_ylabel('Normalized Intensity')
    ax2.set_xlim(300, 3100)
    ax2.set_ylim((14000, 24000))

    ax3.set_title('Flattened, Masked Spectra')
    ax3.set_xlabel('Register Number')
    ax3.set_ylabel('Normalized Intensity')
    ax3.set_xlim(300, 3100)
    ax3.set_ylim((14000, 24000))

    ax4.set_title('Flattened, Masked, Interpolated Spectra')
    ax4.set_xlabel('Register Number')
    ax4.set_ylabel('Normalized Intensity')
    ax4.set_xlim(300, 3100)
    ax4.set_ylim((14000, 24000))

    ax5.set_title('Flattened, Masked, Interpolated, Individually Normalized Spectra')
    ax5.set_xlabel('Register Number')
    ax5.set_ylabel('Normalized Intensity')
    ax5.set_xlim(300, 3100)
    ax5.set_ylim(0.8, 2)

    ax6.set_title('Flattened, Masked, Interpolated, Individually Normalized Spectra')
    ax6.set_xlabel('Register Number')
    ax6.set_ylabel('Normalized Intensity')
    ax6.set_xlim(1500, 2000)
    ax6.set_ylim(0.8, 2)

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()
    fig6.tight_layout()

    name = tab[-8:-4]

    fig1.savefig(project_path + name + "raw_spectra.png")
    fig2.savefig(project_path + name + "flat_spectra.png")
    fig3.savefig(project_path + name + "flat_masked_spectra.png")
    fig4.savefig(project_path + name + "flat_masked_interpolated_spectra.png")
    fig5.savefig(project_path + name + "flat_masked_interpolated_indi_norm_spectra.png")
    fig6.savefig(project_path + name + "Final Spectra range.png")

# Show all plots
plt.show()

# ----- Average spectra in bunches of 109 spectra, corresponding to 5 degrees of sky movement
# average batches of 109 spectra
batch_size = 109  # batch size

# loop for processing the datafiles with the spectra sums.
for tab in filelist:

    with open(sum_path + tab + '_sums.txt', 'r') as sum_tab:

        ra = []
        dec = []
        signal = []
        counter = 0
        for line in sum_tab:

            spec_line = line.split()

            ra.append(float(spec_line[0]))
            dec.append(float(spec_line[1]))
            signal.append(float(spec_line[2]))

            if len(signal) == batch_size:
                # calulate mean and standard deviation for batches of 109 entries
                ra_mean = np.mean(np.array(ra))
                ra_std = np.std(np.array(ra))

                dec_mean = np.mean(np.array(dec))
                dec_std = np.std(np.array(dec))

                signal_mean = np.mean(np.array(signal))
                signal_std = np.std(np.array(signal))

                with open(project_path + tab + '_means.txt', 'a') as mean_tab:
                    mean_tab.write(
                        "%s\t%s\t%s\t%s\t%s\t%s\n" % (ra_mean, ra_std, dec_mean, dec_std, signal_mean, signal_std))
                ra = []
                dec = []
                signal = []

                counter += 1
