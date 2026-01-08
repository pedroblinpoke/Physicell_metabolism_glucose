from scipy.stats import qmc
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from physicool.optimization import PhysiCellBlackBox
import pandas as pd
from scipy import io as sio
from pathlib import Path
from scipy.spatial import ConvexHull
from typing import List, Dict
import math
from sklearn.cluster import DBSCAN
import subprocess


def classify_in_clusters(cells):
    # Combine position arrays into a single 2D array
    ids = cells['ID']
    positions_x = cells['position_x']
    positions_y = cells['position_y']
    positions_z = cells['position_z']
    #print(f"Cells:{cells}")
    positions = np.column_stack((positions_x, positions_y, positions_z))
    # Perform DBSCAN clustering
    clustering = DBSCAN(eps=190, min_samples=2).fit(positions)

    # Convert cells dictionary to DataFrame
    cells_df = pd.DataFrame(cells)

    # Add 'cluster' column to DataFrame
    cells_df['cluster'] = clustering.labels_

    #print(f"Cells DataFrame: {cells_df}")
    return cells_df
def plot_spheres(dataframe, radius=5):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    clusters = dataframe['cluster'].unique()

    for cluster in clusters:
        cluster_data = dataframe[dataframe['cluster'] == cluster]
        x = cluster_data['position_x']
        y = cluster_data['position_y']
        z = cluster_data['position_z']

        ax.scatter(x, y, z, s=4 * math.pi*(radius**3), label=f'Cluster {cluster}')

    ax.set_xlabel('Position X')
    ax.set_ylabel('Position Y')
    ax.set_zlabel('Position Z')
    ax.legend()
    ax.set_title('3D Scatter Plot with Spheres')

    plt.show()


def update_SBML(new_k_gly: float, new_k_aer: float, new_k_ana: float, new_k_ene: float, new_k_prolif: float, new_death_thresh: float, new_prolif_thresh: float, old_params):
    xml_file_path = 'C:/Users/Usuario/Documents/GitHub/PhysiCell_base_velocity/config/demo.xml'
    xml_file_path_write = 'C:/Users/Usuario/Documents/GitHub/PhysiCell_base_velocity/config'

    old_k_glu = old_params['old_k_glu']
    old_k_aer = old_params['old_k_aer']
    old_k_ana = old_params['old_k_ana']
    old_k_ene = old_params['old_k_ene']
    old_atp_sup = old_params['old_atp_sup']
    old_atp_inf = old_params['old_atp_inf']
    old_k_prolif = old_params['old_k_prolif']

    """
    print("k glu:", new_k_glu)
    print("k aer:", new_k_aer)
    print("k ana:", new_k_ana)
    print("k ene:", new_k_ene)
    print("k prolif:", new_k_prolif)
    print("atp sup:", new_atp_sup)
    print("atp inf:", new_atp_inf)
    """

    k_gly = '<parameter id="k_gly" name="k_gly" value="' + str(old_k_glu) + '" constant="true"/>'
    k_aer = '<parameter metaid="COPASI10" id="k_aer" name="k_aer" value="' + str(old_k_aer) + '" constant="true">'
    k_ana = '<parameter metaid="COPASI11" id="k_ane" name="k_ana" value="' + str(old_k_ana) + '" constant="true">'
    k_ene = '<parameter metaid="COPASI12" id="k_usage" name="k_ene" value="' + str(old_k_ene) + '" constant="true">'
    k_prolif = '<parameter id="k_prolif" name="k_prolif" value="' + str(old_k_prolif) + '" constant="true"/>'
    death_thresh = '<parameter metaid="COPASI13" id="energy_death_thresh" name="energy_death_thresh" value="' + str(old_atp_inf) + '" constant="true">'
    prolif_thresh = '<parameter metaid="COPASI14" id="energy_prolif_thresh" name="energy_prolif_thresh" value="' + str(old_atp_sup) + '" constant="true">'

    # Find the position in number of the value attribute
    start_index_k_gly = k_gly.find('value="') + len('value="')
    end_index_k_gly = k_gly.find('" constant', start_index_k_gly)
    start_index_k_aer = k_aer.find('value="') + len('value="')
    end_index_k_aer = k_aer.find('" constant', start_index_k_aer)
    start_index_k_ana = k_ana.find('value="') + len('value="')
    end_index_k_ana = k_ana.find('" constant', start_index_k_ana)
    start_index_k_ene = k_ene.find('value="') + len('value="')
    end_index_k_ene = k_ene.find('" constant', start_index_k_ene)
    start_index_k_prolif = k_prolif.find('value="') + len('value="')
    end_index_k_prolif = k_prolif.find('" constant', start_index_k_prolif)
    start_index_death_thresh = death_thresh.find('value="') + len('value="')
    end_index_death_thresh = death_thresh.find('" constant', start_index_death_thresh)
    start_index_prolif_thresh = prolif_thresh.find('value="') + len('value="')
    end_index_prolif_thresh = prolif_thresh.find('" constant', start_index_prolif_thresh)

    # Open the XML file and read its lines
    with open(xml_file_path, 'r') as file:
        lines = file.readlines()

    # Iterate over the lines and modify variables
    modified_lines = []
    for line in lines:
        if k_gly in line:
            # Replace the old value with the new value
            modified_syntax = k_gly[:start_index_k_gly] + str(new_k_gly) + k_gly[end_index_k_gly:]
            #print("Modified syntax",modified_syntax)
            # Update the value
            modified_lines.append("\t \t \t" + modified_syntax + "\n")
        elif k_aer in line:
            # Replace the old value with the new value
            modified_syntax = k_aer[:start_index_k_aer] +  str(new_k_aer) + k_aer[end_index_k_aer:]
            #print("Modified syntax", modified_syntax)
            # Update the value
            modified_lines.append("\t \t \t" + modified_syntax + "\n")
        elif k_ana in line:
            # Replace the old value with the new value
            modified_syntax = k_ana[:start_index_k_ana] +  str(new_k_ana) + k_ana[end_index_k_ana:]
            #print("Modified syntax", modified_syntax)
            # Update the value
            modified_lines.append("\t \t \t" + modified_syntax + "\n")
        elif k_ene in line:
            # Replace the old value with the new value
            modified_syntax = k_ene[:start_index_k_ene] +  str(new_k_ene) + k_ene[end_index_k_ene:]
            #print("Modified syntax", modified_syntax)
            # Update the value
            modified_lines.append("\t \t \t" + modified_syntax + "\n")
        elif k_prolif in line:
            # Replace the old value with the new value
            modified_syntax = k_prolif[:start_index_k_prolif] +  str(new_k_prolif) + k_prolif[end_index_k_prolif:]
            #print("Modified syntax", modified_syntax)
            #Update the value
            modified_lines.append("\t \t \t" + modified_syntax + "\n")
        elif death_thresh in line:
            # Replace the old value with the new value
            modified_syntax = death_thresh[:start_index_death_thresh] +  str(new_death_thresh) + death_thresh[
                                      end_index_death_thresh:]
            #print("Modified syntax", modified_syntax)
            # Update the value
            modified_lines.append("\t \t \t" + modified_syntax + "\n")
        elif prolif_thresh in line:
            # Replace the old value with the new value
            modified_syntax = prolif_thresh[:start_index_prolif_thresh] +  str(new_prolif_thresh) + prolif_thresh[
                                                                                              end_index_prolif_thresh:]
            #print("Modified syntax", modified_syntax)
            # Update the value
            modified_lines.append("\t \t \t" + modified_syntax + "\n")
        else:
            modified_lines.append(line)

    # Define the file path within the folder

    file_path = os.path.join(xml_file_path_write, 'demo.xml')

    # Write the modified content back to the file
    try:
        with open(file_path, 'w') as file:
            file.writelines(modified_lines)
        #print("XML file successfully created.")
    except Exception as e:
        print("An error occurred while creating the XML file:", str(e))


def shoelace(x_y):
    x_y = np.array(x_y)
    x_y = x_y.reshape(-1,2)

    x = x_y[:,0]
    y = x_y[:,1]

    S1 = np.sum(x*np.roll(y,-1))
    S2 = np.sum(y*np.roll(x,-1))

    area = .5*np.absolute(S1 - S2)

    return area

def get_cell_data(timestep: float, folder_name: Path):
    """Returns a dictionary with the cell output data for the selected variables."""
    #print(f"Timestep:{timestep}")
    #print(f"Folder name:{folder_name}")
    #print(f"Variables:{variables}")
    # All possible output variables written by PhysiCell



    variables = ['ID','current_phase']

    data_labels = [
        'ID',
        'position_x', 'position_y', 'position_z',
        'total_volume',
        'cell_type',
        'cycle_model', 'current_phase', 'elapsed_time_in_phase',
        'nuclear_volume', 'cytoplasmic_volume',
        'fluid_fraction', 'calcified_fraction',
        'orientation_x', 'orientation_y', 'orientation_z',
        'polarity',
        'migration_speed',
        'motility_vector_x', 'motility_vector_y', 'motility_vector_z',
        'migration_bias',
        'motility_bias_direction_x', 'motility_bias_direction_y', 'motility_bias_direction_z',
        'persistence_time',
        'motility_reserved'
    ]

    # Create path name
    time_str = str(timestep).zfill(8)
    file_name = 'output{}_cells_physicell.mat'.format(time_str)
    path_name = folder_name / file_name
    #print(path_name)

    # Read output file
    #print(sio.loadmat(path_name))
    cell_data = sio.loadmat(path_name)['cells']
    data = sio.loadmat(path_name)
    #print(f"Cell data:{cell_data}")
    #print(f"Data:{sio.loadmat(path_name)}")
    # Select and save the variables of interest
    variables_indexes = [data_labels.index(var) for var in variables]
    #print(f"Variable indexes:{variables_indexes}")
    #print(f"Variable :{variables}")
    cells = {var: cell_data[index, :]
             for var, index in zip (variables, variables_indexes)}
    #print(cells)
    return cells



def get_death_cell_percentage_all_time(timemax , folder_name):
    timestep = 0
    # save_every=60
    day_0 = 0
    day_1 = 24
    day_2 = 48
    day_3 = 72
    day_4 = 96
    day_5 = 120
    day_6 = 144
    day_7 = 168

    days_to_save = [day_0, day_1, day_2, day_3, day_4, day_5, day_6, day_7]
    # days_to_save = [day_0, day_1, day_2]
    comp_df = pd.DataFrame(
        columns=['Time', 'Size_Tot', 'Size_not_1_or_2', 'Num_of_1', 'Num_of_2', 'Nums_of_spheroids_and_cells',
                 'Nums_of_spheroids'])
    while timestep < total_time + 1:
        if timestep in days_to_save:
            # print(f"Timestep:{timestep}")
            spheroids_areas = calculate_area_ConvexHull(timestep, folder_name)
            spheroids_areas = [value for value in spheroids_areas.values()]
            spheroids_areas_array = np.array(spheroids_areas)
            if spheroids_areas_array.size == 0:
                spheroids_areas_array = np.array([0])
            # print(f"The areas are:{spheroids_areas_array}")
            spheroids_areas_no_1_2_array = spheroids_areas_array[
                (spheroids_areas_array != 1) & (spheroids_areas_array != 2)]
            if spheroids_areas_no_1_2_array.size == 0:
                spheroids_areas_no_1_2_array = np.array([0])
            # print(f"The areas with no 1 or 2 are:{spheroids_areas_no_1_2_array}")

            comp_df.loc[timestep, 'Time'] = timestep

            if spheroids_areas_array.size == 0:
                comp_df.loc[timestep, 'Size_Tot'] = 0
            if spheroids_areas_array.size == 1:
                comp_df.loc[timestep, 'Size_Tot'] = spheroids_areas_array
            else:
                comp_df.loc[timestep, 'Size_Tot'] = sum(spheroids_areas_array) / len(spheroids_areas_array)

            if spheroids_areas_no_1_2_array.size == 0:
                comp_df.loc[timestep, 'Size_not_1_or_2'] = 0
            if spheroids_areas_no_1_2_array.size == 1:
                comp_df.loc[timestep, 'Size_not_1_or_2'] = spheroids_areas_no_1_2_array
            else:
                comp_df.loc[timestep, 'Size_not_1_or_2'] = sum(spheroids_areas_no_1_2_array) / len(
                    spheroids_areas_no_1_2_array)

            comp_df.loc[timestep, 'Num_of_1'] = count_ones_in_list(spheroids_areas)

            comp_df.loc[timestep, 'Num_of_2'] = count_two_in_list(spheroids_areas)

            comp_df.loc[timestep, 'Nums_of_spheroids_and_cells'] = count_the_clusters_and_cells(timestep, folder_name)

            comp_df.loc[timestep, 'Nums_of_spheroids'] = count_the_clusters(timestep, folder_name)

        timestep = timestep + 1
    # print(f"Compf:{comp_df}")
    return comp_df
def build_path_name(timestep, folder_name, data_type):
    """Returns a Path object with the adequate PhysiCell output name format.

    Uses the standard PhysiCell output filename structure,
    (output{data_type}_{time_id}.mat)

    Parameters
    ----------
    timestep : int
        The time point at which the output was recorded
    folder_path: Path
        The path to the folder where the output (.mat, .xml) files are stored
    data_type: string
        The type of data to be retrieved (cells for cell-based data and
        microenvironment for the continuum variables)
    """

    # All possible file types written by PhysiCell
    data_type_name = {
        'cells': 'cells_physicell'
    }

    # Variable definition
    file_name = data_type_name[data_type]
    time_str = str(timestep).zfill(8)
    file_name = 'output{}_{}.mat'.format(time_str, file_name)

    # Build path name
    path_name = folder_name / file_name

    return path_name

def get_cell_positions(timestep, folder_name):
    """Returns a dictionary with the cell output data for the given variables.

    Parameters
    ----------
    timestep : int
        The time point at which the output was recorded
    folder_path: Path
        The path to the folder where the output (.mat, .xml) files are stored
    variables : list
        The variables to be extracted from the output files. If variables
        are not defined, all the available outputs will be saved.
    """

    variables = ['ID', 'position_x', 'position_y', 'position_z']
    # All possible output variables written by PhysiCell
    data_labels = ['ID','position_x', 'position_y', 'position_z']

    # Variable definition
    cells = {}
    file_type = 'cells'

    # Build path name
    path_name = build_path_name(timestep, folder_name, file_type)
    print(path_name)
    # Read output file
    cell_data = sio.loadmat(path_name)['cells']

    # Select and save the variables of interest
    variables_indexes = [data_labels.index(var) for var in variables]

    for index, var, in zip(variables_indexes, variables):
        cells[var] = cell_data[index, :]

    return cells

def get_cell_data_clasify(timestep: float, folder_name: Path):
    """Returns a dictionary with the cell output data for the selected variables and classify them in clusters."""
    #print(f"Timestep:{timestep}")
    #print(f"Folder name:{folder_name}")
    #print(f"Variables:{variables}")
    # All possible output variables written by PhysiCell



    variables = ['ID','position_x', 'position_y', 'position_z']

    data_labels = [
        'ID',
        'position_x', 'position_y', 'position_z',
        'total_volume',
        'cell_type',
        'cycle_model', 'current_phase', 'elapsed_time_in_phase',
        'nuclear_volume', 'cytoplasmic_volume',
        'fluid_fraction', 'calcified_fraction',
        'orientation_x', 'orientation_y', 'orientation_z',
        'polarity',
        'migration_speed',
        'motility_vector_x', 'motility_vector_y', 'motility_vector_z',
        'migration_bias',
        'motility_bias_direction_x', 'motility_bias_direction_y', 'motility_bias_direction_z',
        'persistence_time',
        'motility_reserved'
    ]

    # Create path name
    time_str = str(timestep).zfill(8)
    file_name = 'output{}_cells_physicell.mat'.format(time_str)
    path_name = folder_name / file_name
    #print(path_name)

    # Read output file
    #print(sio.loadmat(path_name))
    cell_data = sio.loadmat(path_name)['cells']
    data = sio.loadmat(path_name)
    #print(f"Cell data:{cell_data}")
    #print(f"Data:{sio.loadmat(path_name)}")
    # Select and save the variables of interest
    variables_indexes = [data_labels.index(var) for var in variables]
    #print(f"Variable indexes:{variables_indexes}")
    #print(f"Variable :{variables}")
    cells = {var: cell_data[index, :]
             for var, index in zip (variables, variables_indexes)}
    #print(cells)
    return cells
def calculate_area_ConvexHull(timestep, folder_name):
    cells = get_cell_positions(timestep, folder_name)
    cells_positions_df = np.array([cells['position_x'], cells['position_y'], cells['position_z']])
    print(f"Cells positions:{cells_positions_df}")
    area_of_the_spheroid = {}
    if cells_positions_df.size == 0:
        area_of_the_spheroid[0] = 0
    else:
        cells = classify_in_clusters(cells)
        cells_df = pd.DataFrame(cells)
        clustered_cells_df = cells_df.copy()
        singular_cells = cells_df[cells_df['cluster'] == -1]
        cells_in_clusters = cells_df[cells_df['cluster'] != -1]
        print(f"Celulas en cluster:{cells_in_clusters}")
        #print(f"Celulas solas:{singular_cells}")
        #clustered_cells_df['cluster'] = np.nan
        different_clusters = divide_dataframes_by_cluster(cells_in_clusters)
        #print(f"dif_clusters{different_clusters}")
        total_number_clusters = cells_df['cluster'].max() + 1  # Add 1 because range is exclusive on the upper bound
        total_number_ones =(cells_df['cluster'] == -1).sum()
        print(f"Total number of clusters:{total_number_clusters}")

        count=0

        for num_cell in range(total_number_ones):
            area_of_the_spheroid[count] = 1
            count = count + 1
            #print(f"Count:{count} Area:{area_of_the_spheroid}")

        for num_cluster in range(total_number_clusters):
            #print(count)
            if (different_clusters[num_cluster]['cluster'] == -1).any():
                area_of_the_spheroid[count] = 1
                #print(f"Count:{count} Area:{area_of_the_spheroid}")
            else:
                cells_x = different_clusters[num_cluster]['position_x'].astype(float)
                cells_y = different_clusters[num_cluster]['position_y'].astype(float)
                #print(f"Cells_x:{cells_x}")
                #print(f"Cells_y:{cells_y}")
                # Filter out non-numeric values and convert to floats
                cells_x = [float(x) for x in cells_x]
                cells_y = [float(y) for y in cells_y]
                """
                print("Position X: Position Y:")
                for x, y in zip(cells_x, cells_y):
                    print(f"{x} {y}")
                """
                points_of_the_pol = list(zip(cells_x, cells_y))
                if not points_of_the_pol:  # Check if points_of_the_pol is empty
                    print("0")
                    area_of_the_spheroid[count] = 0
                    #print(f"Count:{count} Area:{area_of_the_spheroid}")
                elif len(points_of_the_pol) == 1:
                    print("1")
                    area_of_the_spheroid[count] = 1
                    #print(f"Count:{count} Area:{area_of_the_spheroid}")
                elif len(points_of_the_pol) == 2:
                    print("2")
                    area_of_the_spheroid[count] = 2
                    #print(f"Count:{count} Area:{area_of_the_spheroid}")
                else:
                    print("More than 2")
                    hull = ConvexHull(points_of_the_pol)
                    hull_vertices = [points_of_the_pol[i] for i in hull.vertices]
                    hull_vertices_df = pd.DataFrame(hull_vertices, columns=[ "position_x", "position_y"])
                    print(f"Hull:{hull_vertices_df}")

                    i = np.arange(len(hull_vertices_df))
                    area_of_the_spheroid[count] = shoelace(hull_vertices_df)
                    #print(f"Count:{count} Area:{area_of_the_spheroid}")
                    #print(f"Area shoelace:{area_of_the_spheroid}")


                    #Plot points and perymeter
                    #plt.plot(*zip(*points_of_the_pol), 'ro')
                    #plt.plot(*zip(*hull_vertices), 'bo-')
                    #plt.show()
                    #plt.savefig("shoelace.png")


                    #Calculate the area with the shoelace formula
                    #i = np.arange(len(hull_vertices))
                    #Area = np.abs(np.sum(x[i - 1] * y[i] - x[i] * y[i - 1]) * 0.5)
            count=count+1
    print(f"AREAAA:{area_of_the_spheroid}")
    return area_of_the_spheroid
def count_ones_in_list(lst):
    count = 0
    for item in lst:
        if item == 1:
            count += 1
    return count

def count_two_in_list(lst):
    count = 0
    for item in lst:
        if item == 2:
            count += 1
    return count

def count_the_clusters_and_cells(timestep, folder_name):
    cells_data = get_cell_positions(timestep, folder_name)

    if not cells_data:
        # If cells_data is empty, there are no cells at this timestep
        return 0

    # Convert dictionary values to arrays
    cells_positions = np.array([cells_data['position_x'], cells_data['position_y'], cells_data['position_z']])

    # Check if there are no cells
    if cells_positions.size == 0:
        return 0
    elif cells_positions.shape[1] == 1:
        return 1
    else:
        # Now you can continue your logic using cells_positions
        #print(f"DATA:{cells_data}")
        cells_clust = classify_in_clusters(cells_data)


        count_minus_1 = cells_clust[cells_clust['cluster'] == -1].shape[0]

        max_cluster_value = cells_clust['cluster'].max()
        if max_cluster_value == -1:
            max_cluster_value = 0

        return count_minus_1 + max_cluster_value + 1
def count_the_clusters(timestep, folder_name):
    cells_data = get_cell_positions(timestep, folder_name)

    if not cells_data:
        # If cells_data is empty, there are no cells at this timestep
        return 0

    # Convert dictionary values to arrays
    cells_positions = np.array([cells_data['position_x'], cells_data['position_y'], cells_data['position_z']])

    if cells_positions.size == 0:
        return 0
    elif cells_positions.shape[1] == 1:
        return 1
    else:
        cells_clust = classify_in_clusters(cells_data)
        max_cluster_value = cells_clust['cluster'].max()
        if max_cluster_value == -1:
            max_cluster_value=0
        #print("Maximum value in the 'cluster' column:", max_cluster_value)
    return max_cluster_value+1

def count_ind_cells(timestep, folder_name):
    cells_clust = get_cell_positions(timestep, folder_name)
    cells_clust = classify_in_clusters(cells_clust)

    count_minus_1 = cells_clust[cells_clust['cluster'] == -1].shape[0]
    print("Timestep", timestep)
    print("Number of occurrences of -1 in the 'cluster' column:", count_minus_1)

    return count_minus_1

def get_death_cell_percentage(timestep, folder_name):
    cells = get_cell_data(timestep, folder_name)
    #print(cells)
    cells_life = np.array([cells['current_phase']])

    count_live_cells = np.sum(cells_life == 100)
    total_elements = cells_life.size
    percentage = (count_live_cells / total_elements) * 100
    #print(percentage)
    return percentage
def calculate_all_areas(total_time, folder_name):
    timestep = 0
    #save_every=60
    day_0=0
    day_1=1
    day_2=2
    day_3=3
    day_4=4
    day_5=5
    day_6=6
    day_7=7

    days_to_save=[day_0, day_1, day_2, day_3, day_4, day_5, day_6, day_7]
    #days_to_save = [day_0, day_1, day_2]
    comp_df = pd.DataFrame(columns=['Time', 'Size_Tot', 'Size_not_1_or_2', 'Num_of_1', 'Num_of_2', 'Nums_of_spheroids_and_cells', 'Nums_of_spheroids', 'Percentage_of_death'])
    while timestep < total_time+1:
        if timestep in days_to_save:
            print(f"Timestep:{timestep}")
            spheroids_areas =  calculate_area_ConvexHull(timestep, folder_name)
            spheroids_areas = [value for value in spheroids_areas.values()]
            spheroids_areas_array= np.array(spheroids_areas)
            if spheroids_areas_array.size == 0:
                spheroids_areas_array = np.array([0])
            #print(f"The areas are:{spheroids_areas_array}")
            spheroids_areas_no_1_2_array = spheroids_areas_array[(spheroids_areas_array != 1) & (spheroids_areas_array != 2)]
            if spheroids_areas_no_1_2_array.size == 0:
                spheroids_areas_no_1_2_array = np.array([0])
            #print(f"The areas with no 1 or 2 are:{spheroids_areas_no_1_2_array}")

            comp_df.loc[timestep, 'Time'] = timestep

            if spheroids_areas_array.size == 0:
                comp_df.loc[timestep, 'Size_Tot'] = 0
            if spheroids_areas_array.size == 1:
                comp_df.loc[timestep, 'Size_Tot'] = spheroids_areas_array
            else:
                comp_df.loc[timestep, 'Size_Tot'] = sum(spheroids_areas_array)/len(spheroids_areas_array)

            if spheroids_areas_no_1_2_array.size == 0:
                comp_df.loc[timestep, 'Size_not_1_or_2'] = 0
            if spheroids_areas_no_1_2_array.size == 1:
                comp_df.loc[timestep, 'Size_not_1_or_2'] = spheroids_areas_no_1_2_array
            else:
                comp_df.loc[timestep, 'Size_not_1_or_2'] = sum(spheroids_areas_no_1_2_array)/len(spheroids_areas_no_1_2_array)

            comp_df.loc[timestep, 'Num_of_1'] = count_ones_in_list(spheroids_areas)

            comp_df.loc[timestep, 'Num_of_2'] = count_two_in_list(spheroids_areas)

            comp_df.loc[timestep, 'Nums_of_spheroids_and_cells'] = count_the_clusters_and_cells(timestep, folder_name)

            comp_df.loc[timestep, 'Nums_of_spheroids'] = count_the_clusters(timestep, folder_name)

            comp_df.loc[timestep, 'Percentage_of_death'] = get_death_cell_percentage(timestep, folder_name)

        timestep = timestep+1
    #print(f"Compf:{comp_df}")
    return comp_df
def divide_dataframes_by_cluster(dataframe):
    # Get unique cluster values
    clusters = dataframe['cluster'].unique()

    # Create an empty dictionary to store DataFrames for each cluster
    cluster_dataframes = {}

    # Iterate over unique cluster values
    for cluster in clusters:
        # Filter the original DataFrame for the current cluster
        cluster_data = dataframe[dataframe['cluster'] == cluster]

        # Store the filtered DataFrame in the dictionary
        cluster_dataframes[cluster] = cluster_data

    return cluster_dataframes


new_directory = Path(f'C:/Users/Usuario/Documents/GitHub/PhysiCell_base_velocity')
os.chdir(new_directory)
# 7 params
#Create the parameter configuration
number_of_parameters=7
sampler = qmc.Sobol(d=number_of_parameters, scramble=False)
sample = sampler.random_base2(m=9)
num_of_params_config = sample.size/7
#print(sample)
print(f"Numero de configuraciones:{num_of_params_config}")
#sample

#Separate the list of parameters
k_glu = sample[:, 0]
k_aer = sample[:, 1]
k_ana = sample[:, 2]
k_ene = sample[:, 3]
k_prolif = sample[:, 4]
atp_sup = sample[:, 5]
atp_inf = sample[:, 6]

# Rescale the parameters
k_glu = (1e-11 - 1e-15) * k_glu + 1e-15
k_aer = (1e-11 - 1e-14) * k_aer + 1e-14
k_ana = (1e-11 - 1e-14) * k_ana + 1e-14
k_ene = (1e-1 - 1e-4) * k_ene + 1e-5
k_prolif = ((0.00077) - (0.00045)) * k_prolif + (0.00045)
atp_sup = (1300 - 800) * atp_sup + 800
atp_inf = (800 - 300) * atp_inf + 300



print("k glu:", k_glu)
print("k aer:", k_aer)
print("k ana:", k_ana)
print("k ene:", k_ene)
print("k prolif:", k_prolif)
print("atp sup:", atp_sup)
print("atp inf:", atp_inf)


"""
n=0
m=171
k_glu = k_glu[n:m]
k_aer = k_aer[n:m]
k_ana = k_ana[n:m]
k_ene = k_ene[n:m]
k_prolif = k_prolif[n:m]
atp_sup = atp_sup[n:m]
atp_inf = atp_inf[n:m]

print("k glu:", k_glu)
print("k aer:", k_aer)
print("k ana:", k_ana)
print("k ene:", k_ene)
print("k prolif:", k_prolif)
print("atp sup:", atp_sup)
print("atp inf:", atp_inf)
"""
#plot the parameter configuration
fig, axes = plt.subplots(number_of_parameters, number_of_parameters, figsize=(15, 15))

data = [k_glu, k_aer, k_ana, k_ene, k_prolif, atp_sup, atp_inf]
labels = ['k_glu', 'k_aer', 'k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf']

subprocess.run(["make", "data-cleanup"], check=True)
#Plot the parameters value two by two
"""
for i in range(number_of_parameters):
    for j in range(number_of_parameters):
        if i != j:
            axes[i, j].scatter(data[i], data[j])
            axes[i, j].set_xlabel(labels[i])
            axes[i, j].set_ylabel(labels[j])

# Adjust layout
plt.tight_layout()
plt.show()
"""


##### The parameters on demo.xml has to be set to 0 before the simulation manually
old_params = {'old_k_glu': 0.0,
              'old_k_aer': 0.0,
              'old_k_ana': 0.0,
              'old_k_ene': 0.0,
              'old_atp_sup': 0.0,
              'old_atp_inf': 0.0,
              'old_k_prolif': 0.0}

#Initialice the final dataframe
final_DataFrame_area= pd.DataFrame(columns=['k_glu', 'k_aer','k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf',
                                       'AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7'])
final_DataFrame_area_no12= pd.DataFrame(columns=['k_glu', 'k_aer','k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf',
                                       'AreaD0_no12', 'AreaD1_no12', 'AreaD2_no12', 'AreaD3_no12', 'AreaD4_no12', 'AreaD5_no12', 'AreaD6_no12', 'AreaD7_no12'])
final_DataFrame_ones= pd.DataFrame(columns=['k_glu', 'k_aer','k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf',
                                       'Ones_D0', 'Ones_D1', 'Ones_D2', 'Ones_D3', 'Ones_D4', 'Ones_D5', 'Ones_D6', 'Ones_D7'])
final_DataFrame_twos= pd.DataFrame(columns=['k_glu', 'k_aer','k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf',
                                       'Two_D0', 'Two_D1', 'Two_D2', 'Two_D3', 'Two_D4', 'Two_D5', 'Two_D6', 'Two_D7'])
final_DataFrame_spheroids_cells= pd.DataFrame(columns=['k_glu', 'k_aer','k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf',
                                       'Sp_cells_D0', 'Sp_cells_D1', 'Sp_cells_D2', 'Sp_cells_D3', 'Sp_cells_D4', 'Sp_cells_D5', 'Sp_cells_D6', 'Sp_cells_D7'])
final_DataFrame_spheroids= pd.DataFrame(columns=['k_glu', 'k_aer','k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf',
                                       'Sp_D0', 'Sp_D1', 'Sp_D2', 'Sp_D3', 'Sp_D4', 'Sp_D5', 'Sp_D6', 'Sp_D7'])
final_DataFrame_deaths= pd.DataFrame(columns=['k_glu', 'k_aer','k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf',
                                       'Death_D0', 'Death_D1', 'Death_D2', 'Death_D3', 'Death_D4', 'Death_D5', 'Death_D6', 'Death_D7'])

#For that goes through all the parameter configurations
i=0

output_files = {
    'areas': 'output_areas.txt',
    'areas_no12': 'output_areas_no12.txt',
    'ones': 'output_ones.txt',
    'twos': 'output_twos.txt',
    'spheroids_cells': 'output_spheroids_cells.txt',
    'spheroids': 'output_spheroids.txt',
    'deaths': 'output_deaths.txt'
}

for file in output_files.values():
    with open(file, 'w') as f:
        pass
#num_of_params_config=172
for i in range(int(num_of_params_config)):
    print("Parameter configuration number:",i)
    #Store the params in variables

    new_k_glu = float(k_glu[i])
    new_k_aer = float(k_aer[i])
    new_k_ana = float(k_ana[i])
    new_k_ene = float(k_ene[i])
    new_k_prolif = float(k_prolif[i])
    new_atp_sup = float(atp_sup[i])
    new_atp_inf = float(atp_inf[i])


    #Print the parameters

    print("k glu:", new_k_glu)
    print("k aer:", new_k_aer)
    print("k ana:", new_k_ana)
    print("k ene:", new_k_ene)
    print("k prolif:", new_k_prolif)
    print("atp sup:", new_atp_sup)
    print("atp inf:", new_atp_inf)


    update_SBML(new_k_glu, new_k_aer, new_k_ana, new_k_ene, new_k_prolif, new_atp_inf, new_atp_sup, old_params)

    old_params = {'old_k_glu': new_k_glu,
                  'old_k_aer': new_k_aer,
                  'old_k_ana': new_k_ana,
                  'old_k_ene': new_k_ene,
                  'old_k_prolif': new_k_prolif,
                  'old_atp_inf': new_atp_inf,
                  'old_atp_sup': new_atp_sup}



    num_simulations = 1  # Number of simulations to process
    print(f"Pre Data-cleanup")

    os.system("make data-cleanup")
    print(f"Pre Make")

    os.system("make")
    print(f"Empieza simulacion:{i}")
    my_model = PhysiCellBlackBox(project_name="ode_energy")
    my_model.run(number_of_replicates=num_simulations, keep_files=True)
    print(f"Termina simulacion:{i}")
    total_time_simulations =7 #Tiempo total en pasos temporales
    areas_data = {}  # Dictionary to store areas data for each simulation
    folder_name = Path(f'C:/Users/Usuario/Documents/GitHub/PhysiCell_base_velocity/temp')

    #print(folder_name)
    current_areas_data=calculate_all_areas(total_time_simulations, folder_name)

    #print(f"Areas_data:{current_areas_data}")

    df_to_append_areas_tot = {
        'k_glu': new_k_glu, 'k_aer': new_k_aer, 'k_ana': new_k_ana, 'k_ene': new_k_ene, 'k_prolif': new_k_prolif,
        'atp_sup': new_atp_sup, 'atp_inf': new_atp_inf, 'atp_sup-atp_inf': new_atp_sup - new_atp_inf,
        'AreaD0': current_areas_data.loc[0, 'Size_Tot'], 'AreaD1': current_areas_data.loc[1, 'Size_Tot'],
        'AreaD2': current_areas_data.loc[2, 'Size_Tot'], 'AreaD3': current_areas_data.loc[3, 'Size_Tot'],
        'AreaD4': current_areas_data.loc[4, 'Size_Tot'], 'AreaD5': current_areas_data.loc[5, 'Size_Tot'],
        'AreaD6': current_areas_data.loc[6, 'Size_Tot'], 'AreaD7': current_areas_data.loc[7, 'Size_Tot']
    }

    df_to_append_areas_no12 = {
        'k_glu': new_k_glu, 'k_aer': new_k_aer, 'k_ana': new_k_ana, 'k_ene': new_k_ene, 'k_prolif': new_k_prolif,
        'atp_sup': new_atp_sup, 'atp_inf': new_atp_inf, 'atp_sup-atp_inf': new_atp_sup - new_atp_inf,
        'AreaD0no12': current_areas_data.loc[0, 'Size_not_1_or_2'],
        'AreaD1no12': current_areas_data.loc[1, 'Size_not_1_or_2'],
        'AreaD2no12': current_areas_data.loc[2, 'Size_not_1_or_2'],
        'AreaD3no12': current_areas_data.loc[3, 'Size_not_1_or_2'],
        'AreaD4no12': current_areas_data.loc[4, 'Size_not_1_or_2'],
        'AreaD5no12': current_areas_data.loc[5, 'Size_not_1_or_2'],
        'AreaD6no12': current_areas_data.loc[6, 'Size_not_1_or_2'],
        'AreaD7no12': current_areas_data.loc[7, 'Size_not_1_or_2']
    }

    df_to_append_ones = {
        'k_glu': new_k_glu, 'k_aer': new_k_aer, 'k_ana': new_k_ana, 'k_ene': new_k_ene, 'k_prolif': new_k_prolif,
        'atp_sup': new_atp_sup, 'atp_inf': new_atp_inf, 'atp_sup-atp_inf': new_atp_sup - new_atp_inf,
        'Numer_One_D0': current_areas_data.loc[0, 'Num_of_1'], 'Numer_One_D1': current_areas_data.loc[1, 'Num_of_1'],
        'Numer_One_D2': current_areas_data.loc[2, 'Num_of_1'], 'Numer_One_D3': current_areas_data.loc[3, 'Num_of_1'],
        'Numer_One_D4': current_areas_data.loc[4, 'Num_of_1'], 'Numer_One_D5': current_areas_data.loc[5, 'Num_of_1'],
        'Numer_One_D6': current_areas_data.loc[6, 'Num_of_1'], 'Numer_One_D7': current_areas_data.loc[7, 'Num_of_1']
    }

    df_to_append_twos = {
        'k_glu': new_k_glu, 'k_aer': new_k_aer, 'k_ana': new_k_ana, 'k_ene': new_k_ene, 'k_prolif': new_k_prolif,
        'atp_sup': new_atp_sup, 'atp_inf': new_atp_inf, 'atp_sup-atp_inf': new_atp_sup - new_atp_inf,
        'Numer_Two_D0': current_areas_data.loc[0, 'Num_of_2'], 'Numer_Two_D1': current_areas_data.loc[1, 'Num_of_2'],
        'Numer_Two_D2': current_areas_data.loc[2, 'Num_of_2'], 'Numer_Two_D3': current_areas_data.loc[3, 'Num_of_2'],
        'Numer_Two_D4': current_areas_data.loc[4, 'Num_of_2'], 'Numer_Two_D5': current_areas_data.loc[5, 'Num_of_2'],
        'Numer_Two_D6': current_areas_data.loc[6, 'Num_of_2'], 'Numer_Two_D7': current_areas_data.loc[7, 'Num_of_2']
    }

    df_to_append_spheroids_cells = {
        'k_glu': new_k_glu, 'k_aer': new_k_aer, 'k_ana': new_k_ana, 'k_ene': new_k_ene, 'k_prolif': new_k_prolif,
        'atp_sup': new_atp_sup, 'atp_inf': new_atp_inf, 'atp_sup-atp_inf': new_atp_sup - new_atp_inf,
        'Number_spheroids_cells_D0': current_areas_data.loc[0, 'Nums_of_spheroids'],
        'Number_spheroids_cells_D1': current_areas_data.loc[1, 'Nums_of_spheroids'],
        'Number_spheroids_cells_D2': current_areas_data.loc[2, 'Nums_of_spheroids'],
        'Number_spheroids_cells_D3': current_areas_data.loc[3, 'Nums_of_spheroids'],
        'Number_spheroids_cells_D4': current_areas_data.loc[4, 'Nums_of_spheroids'],
        'Number_spheroids_cells_D5': current_areas_data.loc[5, 'Nums_of_spheroids'],
        'Number_spheroids_cells_D6': current_areas_data.loc[6, 'Nums_of_spheroids'],
        'Number_spheroids_cells_D7': current_areas_data.loc[7, 'Nums_of_spheroids']
    }

    df_to_append_spheroids = {
        'k_glu': new_k_glu, 'k_aer': new_k_aer, 'k_ana': new_k_ana, 'k_ene': new_k_ene, 'k_prolif': new_k_prolif,
        'atp_sup': new_atp_sup, 'atp_inf': new_atp_inf, 'atp_sup-atp_inf': new_atp_sup - new_atp_inf,
        'Number_spheroids_D0': current_areas_data.loc[0, 'Nums_of_spheroids'],
        'Number_spheroids_D1': current_areas_data.loc[1, 'Nums_of_spheroids'],
        'Number_spheroids_D2': current_areas_data.loc[2, 'Nums_of_spheroids'],
        'Number_spheroids_D3': current_areas_data.loc[3, 'Nums_of_spheroids'],
        'Number_spheroids_D4': current_areas_data.loc[4, 'Nums_of_spheroids'],
        'Number_spheroids_D5': current_areas_data.loc[5, 'Nums_of_spheroids'],
        'Number_spheroids_D6': current_areas_data.loc[6, 'Nums_of_spheroids'],
        'Number_spheroids_D7': current_areas_data.loc[7, 'Nums_of_spheroids']
    }

    df_to_append_deaths = {
        'k_glu': new_k_glu, 'k_aer': new_k_aer, 'k_ana': new_k_ana, 'k_ene': new_k_ene, 'k_prolif': new_k_prolif,
        'atp_sup': new_atp_sup, 'atp_inf': new_atp_inf, 'atp_sup-atp_inf': new_atp_sup - new_atp_inf,
        'Deaths_D0': current_areas_data.loc[0, 'Percentage_of_death'],
        'Deaths_D1': current_areas_data.loc[1, 'Percentage_of_death'],
        'Deaths_D2': current_areas_data.loc[2, 'Percentage_of_death'],
        'Deaths_D3': current_areas_data.loc[3, 'Percentage_of_death'],
        'Deaths_D4': current_areas_data.loc[4, 'Percentage_of_death'],
        'Deaths_D5': current_areas_data.loc[5, 'Percentage_of_death'],
        'Deaths_D6': current_areas_data.loc[6, 'Percentage_of_death'],
        'Deaths_D7': current_areas_data.loc[7, 'Percentage_of_death']
    }

    # Write each DataFrame to its respective text file
    data_to_write = {
        'areas': df_to_append_areas_tot,
        'areas_no12': df_to_append_areas_no12,
        'ones': df_to_append_ones,
        'twos': df_to_append_twos,
        'spheroids_cells': df_to_append_spheroids_cells,
        'spheroids': df_to_append_spheroids,
        'deaths': df_to_append_deaths
    }

    for key, df_data in data_to_write.items():
        with open(output_files[key], 'a') as file:
            df = pd.DataFrame([df_data])
            df.to_string(file, index=False, header=(i == 0))  # Write header only for the first DataFrame
            file.write('\n')  # Add newline for separation
    os.system('make data-cleanup')

"""
#Tratamiento de datos para su mejor representacion
final_DataFrame_area[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']] = final_DataFrame_area[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']].apply(pd.to_numeric, errors='coerce')
final_DataFrame_area.fillna(0, inplace=True)
final_DataFrame_area = final_DataFrame_area.astype(str).replace('\.0', '', regex=True)
final_DataFrame_area.to_csv('output_areas_modified.txt', index=False, sep='\t')

final_DataFrame_area_no12[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']] = final_DataFrame_area[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']].apply(pd.to_numeric, errors='coerce')
final_DataFrame_area_no12.fillna(0, inplace=True)
final_DataFrame_area_no12 = final_DataFrame_area_no12.astype(str).replace('\.0', '', regex=True)
final_DataFrame_area_no12.to_csv('output_areas_no12_modified.txt', index=False, sep='\t')

final_DataFrame_ones[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']] = final_DataFrame_area[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']].apply(pd.to_numeric, errors='coerce')
final_DataFrame_ones.fillna(0, inplace=True)
final_DataFrame_ones = final_DataFrame_ones.astype(str).replace('\.0', '', regex=True)
final_DataFrame_ones.to_csv('output_ones_modified.txt', index=False, sep='\t')

final_DataFrame_twos[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']] = final_DataFrame_area[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']].apply(pd.to_numeric, errors='coerce')
final_DataFrame_twos.fillna(0, inplace=True)
final_DataFrame_twos = final_DataFrame_twos.astype(str).replace('\.0', '', regex=True)
final_DataFrame_twos.to_csv('output_twos_modified.txt', index=False, sep='\t')

final_DataFrame_spheroids[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']] = final_DataFrame_area[['AreaD0', 'AreaD1', 'AreaD2', 'AreaD3', 'AreaD4', 'AreaD5', 'AreaD6', 'AreaD7']].apply(pd.to_numeric, errors='coerce')
final_DataFrame_spheroids.fillna(0, inplace=True)
final_DataFrame_spheroids = final_DataFrame_spheroids.astype(str).replace('\.0', '', regex=True)
final_DataFrame_spheroids.to_csv('output_spheroids_modified.txt', index=False, sep='\t')


#Write in txt files
#(f"Final_Dataframe:{final_DataFrame}")
final_DataFrame_area.to_csv('output_areas.txt', index=False, sep='\t')
final_DataFrame_area_no12.to_csv('output_areas_no12.txt', index=False, sep='\t')
final_DataFrame_ones.to_csv('output_ones.txt', index=False, sep='\t')
final_DataFrame_twos.to_csv('output_twos.txt', index=False, sep='\t')
final_DataFrame_spheroids.to_csv('output_spheroids.txt', index=False, sep='\t')
final_DataFrame_deaths.to_csv('output_deaths.txt', index=False, sep='\t')
"""


