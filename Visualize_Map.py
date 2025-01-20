import os
import sys
from matplotlib import cm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from numpy import number


def custom_round(number, ndigits=0):
    """Always round off"""
    exp = number * 10 ** ndigits
    if abs(exp) - abs(math.floor(exp)) < 0.5:
        return type(number)(math.floor(exp) / 10 ** ndigits)
    return type(number)(math.ceil(exp) / 10 ** ndigits)

def float_to_ra(ra):
    hours = int(ra)
    minutes = int((ra - hours) * 60)
    #seconds = (ra - hours - minutes / 60) * 3600
    return f"{hours}h {minutes}m" #{seconds:.2f}s"

def read_table(table, path):

    print('Start reading', table)
    map_tab = path + '/' + table
    data_list = []

    with open(map_tab, 'r') as file:
        for line in file:
            tab_line = line.split()

            #if only 3 columns are given, change it into 6 columns to be compatible with the other code
            if len(tab_line) == 3:
                if len(tab_line[2]) >= 4:
                    tab_line[2] = tab_line[2][1:4] #if signal has more than 3 characters, remove the first (some values were '<xx.x')

                data_list.append(np.array([tab_line[0], 0, tab_line[1], 0, tab_line[2], 0], dtype=float))

            else:
                data_list.append(np.array(tab_line, dtype=float))

        map_data = np.array(data_list)

        RA, DEC, SIGNAL = map_data[:, 0], map_data[:, 2], map_data[:, 4]

    print('Finished reading', table)
    return RA, DEC, SIGNAL

def fill_matrix(ra_arr, dec_arr, sig_arr, num_file):

    for i in range(len(ra_arr)):

        #if ra larger than 25, print an error
        if ra_arr[i] > 24:
            print('Error:', ra_arr[i])

        #Calculate dec as mean of the input list
        dec = np.mean(dec_arr)
        #Put dec as the first element in the row
        map_mat[num_file][0] = dec

        #Create index for right ascension (RA*10) +1, shift by one to leave the first element as zero
        new_idx = int(custom_round(ra_arr[i] * 10)) + 1

        #if ra >= 24, shift it back to 0
        if new_idx > 239:
            new_idx-=239

        #Put value of right ascension into first row, if not zero, take mean
        if map_mat[0][new_idx] == 0:
            map_mat[0][new_idx] = ra_arr[i]
        else:
            map_mat[0][new_idx] = np.mean(map_mat[0][new_idx]+ra_arr[i])

        # if ra larger than 25, print an error
        if map_mat[0][new_idx] > 25:
            print('Error:', map_mat[0][new_idx])

        map_mat[num_file][new_idx] = sig_arr[i]

    return map_mat

def interpolate_matrix(input_matrix):

    matrix = np.copy(input_matrix)

    for num_row in range(len(matrix)):

        for i in range(1, 240):
            #Do not interpolate if RA is larger than 23.9, i.e. idx > 239
            if i > 239:
                break
            else:

                zero_check = matrix[num_row][i]

                # if first value is 0 leave it to be 0
                if zero_check == 0 and i == 1:
                    matrix[num_row][i] = 0

                # if last value is 0 set it to the previous value
                elif zero_check == 0 and i == len(matrix[num_row]) - 1:
                    matrix[num_row][i] = matrix[num_row][i-1]

                # for all other zeros, set to mean of the two neighbouring non-zero
                elif zero_check == 0 and 0 < i < len(matrix[num_row]) - 1:
                    next_idx = i + 1 #create next idx to current element
                    if next_idx >= len(matrix[num_row]) - 1:
                        break
                    else:
                        while next_idx < len(matrix[num_row]) - 1 and matrix[num_row][next_idx] == 0:
                            next_idx += 1 #increase index as long as value is 0

                        matrix[num_row][i] = np.mean([matrix[num_row][i - 1], matrix[num_row][next_idx]])

        #round elements in the first row (ra)
        for i in range(len(matrix[0])):
            matrix[0][i] = custom_round(matrix[0][i],2)

        #round elements in the first column (dec)
        for i in range(len(matrix)):
            matrix[i][0] = custom_round(matrix[i][0],2)

        print('Row ', num_row, 'done!')
    return matrix


def sort_and_label_matrix(input_matrix):

    matrix = np.copy(input_matrix)

    #create x positions and labels based on position (index) divided by 10
    x_label_pos = np.arange(len(matrix[0]) - 1)
    x_labels = []
    for i in x_label_pos:
        x_labels.append(float_to_ra(int(i)/10))


    #delete first row (right ascension values)
    matrix = np.delete(matrix, 0, axis=0)

    #sort matrix after first column (declination) from highest to lowest value
    mat = matrix[np.argsort(matrix[:, 0])[::-1]]

    #print('idx', np.where(mat[:,0] == -12), np.where(mat[:,0] == 43)) #find positions of declination thresholds from first plot

    #create y ticks and labels
    y_label_pos = np.arange(len(matrix))
    y_labels = mat[:,0]

    #delete first column (declination values)
    final_matrix = np.delete(mat, 0, axis=1)

    return final_matrix, x_label_pos, x_labels, y_label_pos, y_labels


n_file = 1 #counter for files, start at 1

#file path
map_path = "C:/Users/Celina/Masterarbeit/Kartenmessung/First_Final_Data/Map_Data"
map_filelist = os.listdir(map_path)


number_files = len(map_filelist)
map_mat = np.zeros((number_files+1, 241))

for tab in map_filelist:

    #read table and extract ra, dec and signal
    ra_map, dec_map, signal_map = read_table(tab, map_path)

    #fill matrix with ra, dec, and signal
    raw_matrix = fill_matrix(ra_map, dec_map, signal_map, n_file)

    n_file+=1

print('Matrix created')

#interpolate the 'empty' spots
interpolated = interpolate_matrix(raw_matrix)

#sort matrix and delete first row and column after extracting the tick labels
map_matrix, ra_label_pos, ra_labels, dec_label_pos, dec_labels = sort_and_label_matrix(interpolated)

#Create plot
fig, ax = plt.subplots(figsize=(8, 2))

cax = ax.imshow(map_matrix, aspect=4, cmap='plasma')
#fig.colorbar(cax)

#Create tick labels
plt.xticks([ra_label_pos[0], ra_label_pos[30], ra_label_pos[60], ra_label_pos[90], ra_label_pos[120], ra_label_pos[150], ra_label_pos[180], ra_label_pos[210], ra_label_pos[-1]], [ra_labels[0], ra_labels[30], ra_labels[60], ra_labels[90], ra_labels[120], ra_labels[150], ra_labels[180], ra_labels[210], ra_labels[-1]])
plt.yticks([dec_label_pos[0], dec_label_pos[3], dec_label_pos[6], dec_label_pos[9], dec_label_pos[-1]], [dec_labels[0], dec_labels[3], dec_labels[6], dec_labels[9], dec_labels[-1]])

#Create axis labels
plt.xlabel('Right Ascension')
plt.ylabel('Declination [°]')

plt.tight_layout()
plt.savefig("C:/Users/Celina/Masterarbeit/Kartenmessung/Radio_Map_small_area.png")
plt.show()


#----- Cold spot map ------

file_path = "C:/Users/Celina/Masterarbeit/Kartenmessung/"
file = "NVSSCatalog_RA_DEC_FLUX.tab"

#read table
ra_cold, dec_cold, sig_cold = read_table(file, file_path)

#create empty matrix for cold sky map
cold_sky = np.zeros((133, 241))
for j in range(len(ra_cold)):

    #use ra as index for rows
    ra_idx = int(custom_round(ra_cold[j]*10, 0))+1

    #after 24h, start back at index 1
    if ra_idx > 239:
        ra_idx -= 239

    #shift index of declination to 0 (measurement started at -41 deg)
    dec_idx = int(custom_round(dec_cold[j], 0))+41+1

    #fill 0th row with right ascension values
    cold_sky[0][ra_idx] = custom_round(float(ra_cold[j]))

    #fill 0th column with declination values
    cold_sky[dec_idx][0] = custom_round(float(dec_cold[j]))

    #sum all values with the same rounded coordinates
    cold_sky[int(dec_idx)][int(ra_idx)] = (sig_cold[j] + cold_sky[int(dec_idx)][int(ra_idx)])




#sort matrix and delete ra and dec row and column after creating tick labels
cold_mat, ra_cold_pos, ra_cold_label, dec_cold_pos, dec_cold_label= sort_and_label_matrix(cold_sky)

#----- Create figure, cold sky -----
fig, ax = plt.subplots()
cax1 = ax.imshow(cold_mat, cmap='plasma', vmin=0, vmax=8000)

#plot lines to display area of the first plot
ax.axhline(47, color='white', linewidth=0.5, linestyle='--')
ax.axhline(102, color='white', linewidth=0.5, linestyle='--')

#create tick labels
plt.xticks([ra_cold_pos[0],ra_cold_pos[30], ra_cold_pos[60], ra_cold_pos[90], ra_cold_pos[120], ra_cold_pos[150],
            ra_cold_pos[180], ra_cold_pos[210], ra_cold_pos[-1]],
           [ra_cold_label[0], ra_cold_label[30], ra_cold_label[60], ra_cold_label[90], ra_cold_label[120], ra_cold_label[150], ra_cold_label[180], ra_cold_label[210], ra_cold_label[-1]])
plt.yticks(np.array(dec_cold_pos[0::10], dtype=int), np.array(dec_cold_label[0::10], dtype = str))

#create axes label
plt.xlabel('Right Ascension')
plt.ylabel('Declination [°]')

plt.tight_layout()
plt.savefig("C:/Users/Celina/Masterarbeit/Cold_sky.png")
#fig.colorbar(cax1)
plt.show()


