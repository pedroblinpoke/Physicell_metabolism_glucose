import pandas as pd
import re
from pathlib import Path
import plotly.graph_objects as go


input_file = "df_area_no12_mean.txt"
output_file = "df_area_no12_mean_mod.txt"
# Read the input text file
with open(input_file, "r") as file:
    lines = file.readlines()

# Modify the lines to remove square brackets from numbers
modified_lines = []
for line in lines:
    # Find all numbers within square brackets
    numbers_in_brackets = re.findall(r'\[(.*?)\]', line)
    # Replace each number within square brackets with the number itself
    for number_in_brackets in numbers_in_brackets:
        line = line.replace(f"[{number_in_brackets}]", number_in_brackets)
    modified_lines.append(line)

# Write the modified data into a new text file
with open(output_file, "w") as file:
    file.writelines(modified_lines)

print("Modification complete. Output written to", output_file)


#file_path =  Path(f'C:/Users/Usuario/Documents/GitHub/PhysiCell_base_velocity/output_spheroids_modified.txt')


#fig = px.parallel_coordinates(df, color="AreaD1", labels={"AreaD1": "AreaD1",
#                "k_glu": "k_glu", "k_aer": "k_aer",
#                "k_ana": "k_ana", "k_ene": "k_ene", "k_prolif": "k_prolif","atp_sup": "atp_sup","atp_inf": "atp_inf","AreaD1": "AreaD1", "AreaD2": "AreaD2",},
#                             color_continuous_scale=px.colors.diverging.Tealrose,
#                             color_continuous_midpoint=2)
#fig.show()

"""
fig = go.Figure(data=
    go.Parcoords(
        line_color='blue',
        dimensions = list([
            dict(range = [1,5],
                 constraintrange = [1,2], # change this range by dragging the pink line
                 label = 'A', values = [1,4]),
            dict(range = [1.5,5],
                 tickvals = [1.5,3,4.5],
                 label = 'B', values = [3,1.5]),
            dict(range = [1,5],
                 tickvals = [1,2,4,5],
                 label = 'C', values = [2,4],
                 ticktext = ['text 1', 'text 2', 'text 3', 'text 4']),
            dict(range = [1,5],
                 label = 'D', values = [4,2])
        ])
    )
)
fig.show()
"""

file_path_area =  Path(f'C:/Users/Usuario/Documents/GitHub/PhysiCell_base_velocity/df_area_no12_mean.txt')

df_area = pd.read_csv(file_path_area, sep='\t')

print(df_area)

with open(file_path_area, "r") as file:
    first_line = file.readline().strip()
    column_names = first_line.split()
# Step 2: Load the DataFrame using the extracted column names
df_area = pd.read_csv(file_path_area, delim_whitespace=True, skiprows=1, names=column_names)
print("Column names:", column_names[0])

custom_range_AreaD3 = [0, 100000]
custom_range_AreaD5 = [0, 100000]
custom_range_AreaD7 = [0, 100000]

print(df_area)
fig = go.Figure(data=
go.Parcoords(
    line_color='blue',
    dimensions = list([
        dict(range = [float(min(df_area[column_names[0]])), float(max(df_area[column_names[0]]))],
             constraintrange = [float(min(df_area[column_names[0]])), float(max(df_area[column_names[0]]))],
             label = column_names[0], values = df_area[column_names[0]]),
        dict(range = [float(min(df_area[column_names[1]])), float(max(df_area[column_names[1]]))],
             tickvals = [float(min(df_area[column_names[1]])), float(max(df_area[column_names[1]]))],
             label = column_names[1], values = df_area[column_names[1]]),
        dict(range = [float(min(df_area[column_names[2]])), float(max(df_area[column_names[2]]))],
             tickvals = sorted(df_area[column_names[2]].unique()),
             label = column_names[2], values = df_area[column_names[2]]),
        dict(range = [float(min(df_area[column_names[3]])), float(max(df_area[column_names[3]]))],
             label = column_names[3], values = df_area[column_names[3]]),
        dict(range = [float(min(df_area[column_names[4]])), float(max(df_area[column_names[4]]))],
             label = column_names[4], values = df_area[column_names[4]]),
        dict(range = [float(min(df_area[column_names[5]])), float(max(df_area[column_names[5]]))],
             label = column_names[5], values = df_area[column_names[5]]),
        dict(range = [float(min(df_area[column_names[6]])), float(max(df_area[column_names[6]]))],
             label = column_names[6], values = df_area[column_names[6]]),

        dict(range=custom_range_AreaD3,
             label='AreaD3no12', values=df_area['AreaD3no12']),
        dict(range=custom_range_AreaD5,
             label='AreaD5no12', values=df_area['AreaD5no12']),
        dict(range=custom_range_AreaD7,
             label='AreaD7no12', values=df_area['AreaD7no12'])
    ])
)
)
#fig.show()

# Save the plot as an HTML file
fig.write_html('AAAAAA.html')


























input_file = "output_areas.txt"
output_file = "output_areas_modified.txt"
# Read the input text file
with open(input_file, "r") as file:
    lines = file.readlines()

# Modify the lines to remove square brackets from numbers
modified_lines = []
for line in lines:
    # Find all numbers within square brackets
    numbers_in_brackets = re.findall(r'\[(.*?)\]', line)
    # Replace each number within square brackets with the number itself
    for number_in_brackets in numbers_in_brackets:
        line = line.replace(f"[{number_in_brackets}]", number_in_brackets)
    modified_lines.append(line)

# Write the modified data into a new text file
with open(output_file, "w") as file:
    file.writelines(modified_lines)



file_path_area =  Path(f'D:/Resultados_simulaciones/18_06_2024_Dos_cond/Col 25/output_areas.txt')
file_path_area_no12 =  Path(f'D:/Resultados_simulaciones/18_06_2024_Dos_cond/Col 25/output_areas_no12.txt')
file_path_ones =  Path(f'D:/Resultados_simulaciones/18_06_2024_Dos_cond/Col 25/output_ones.txt')
file_path_twos =  Path(f'D:/Resultados_simulaciones/18_06_2024_Dos_cond/Col 25/output_twos.txt')
file_path_spheroids =  Path(f'D:/Resultados_simulaciones/18_06_2024_Dos_cond/Col 25/output_spheroids.txt')

df_area = pd.read_csv(file_path_area, sep='\t')
df_area_no12 = pd.read_csv(file_path_area_no12, sep='\t')
df_ones = pd.read_csv(file_path_ones, sep='\t')
df_twos = pd.read_csv(file_path_twos, sep='\t')
df_spheroids = pd.read_csv(file_path_spheroids, sep='\t')

print(df_area_no12[column_names[0]])

with open(output_file, "r") as file:
    first_line = file.readline().strip()
    column_names = first_line.split()
# Step 2: Load the DataFrame using the extracted column names
df_area = pd.read_csv(output_file, delim_whitespace=True, skiprows=1, names=column_names)
print("Column names:", column_names[0])

# Save the plot as an HTML file
fig.write_html('parallel_plot_area_no12.html')

custom_range_AreaD3 = [0, 18000]
custom_range_AreaD5 = [0, 18000]
custom_range_AreaD7 = [0, 18000]

fig = go.Figure(data=
go.Parcoords(
    line_color='blue',
    dimensions = list([
        dict(range = [float(min(df_area_no12[column_names[0]])), float(max(df_area_no12[column_names[0]]))],
             constraintrange = [float(min(df_area_no12[column_names[0]])), float(max(df_area_no12[column_names[0]]))],
             label = column_names[0], values = df_area_no12[column_names[0]]),
        dict(range = [float(min(df_area_no12[column_names[1]])), float(max(df_area_no12[column_names[1]]))],
             tickvals = [float(min(df_area_no12[column_names[1]])), float(max(df_area_no12[column_names[1]]))],
             label = column_names[1], values = df_area_no12[column_names[1]]),
        dict(range = [float(min(df_area_no12[column_names[2]])), float(max(df_area_no12[column_names[2]]))],
             tickvals = sorted(df_area_no12[column_names[2]].unique()),
             label = column_names[2], values = df_area_no12[column_names[2]]),
        dict(range = [float(min(df_area_no12[column_names[3]])), float(max(df_area_no12[column_names[3]]))],
             label = column_names[3], values = df_area_no12[column_names[3]]),
        dict(range = [float(min(df_area_no12[column_names[4]])), float(max(df_area_no12[column_names[4]]))],
             label = column_names[4], values = df_area_no12[column_names[4]]),
        dict(range = [float(min(df_area_no12[column_names[5]])), float(max(df_area_no12[column_names[5]]))],
             label = column_names[5], values = df_area_no12[column_names[5]]),
        dict(range = [float(min(df_area_no12[column_names[6]])), float(max(df_area_no12[column_names[6]]))],
             label = column_names[6], values = df_area_no12[column_names[6]]),
        dict(range=custom_range_AreaD3,
             label='AreaD3no12', values=df_area_no12['AreaD3no12']),
        dict(range=custom_range_AreaD5,
             label='AreaD5no12', values=df_area_no12['AreaD5no12']),
        dict(range=custom_range_AreaD7,
             label='AreaD7no12', values=df_area_no12['AreaD7no12'])
    ])
)
)
#fig.show()

# Save the plot as an HTML file
fig.write_html('parallel_plot_area_no12.html')

custom_range_AreaD3 = [0, 45]
custom_range_AreaD5 = [0, 45]
custom_range_AreaD7 = [0, 45]

fig = go.Figure(data=
go.Parcoords(
    line_color='blue',
    dimensions = list([
        dict(range = [float(min(df_ones['k_glu'])), float(max(df_ones['k_glu']))],
             constraintrange = [float(min(df_ones['k_glu'])), float(max(df_ones['k_glu']))],
             label = 'k_glu', values = df_ones['k_glu']),
        dict(range = [float(min(df_ones['k_aer'])), float(max(df_ones['k_aer']))],
             tickvals = [float(min(df_ones['k_aer'])), float(max(df_ones['k_aer']))],
             label = 'k_aer', values = df_ones['k_aer']),
        dict(range = [float(min(df_ones['k_ana'])), float(max(df_ones['k_ana']))],
             tickvals = sorted(df_ones['k_ana'].unique()),
             label = 'k_ana', values = df_ones['k_ana']),
        dict(range = [float(min(df_ones['k_ene'])), float(max(df_ones['k_ene']))],
             label = 'k_ene', values = df_ones['k_ene']),
        dict(range = [float(min(df_ones['k_prolif'])), float(max(df_ones['k_prolif']))],
             label = 'k_prolif', values = df_ones['k_prolif']),
        dict(range = [float(min(df_ones['atp_sup'])), float(max(df_ones['atp_sup']))],
             label = 'atp_sup', values = df_ones['atp_sup']),
        dict(range = [float(min(df_ones['atp_inf'])), float(max(df_ones['atp_inf']))],
             label = 'atp_inf', values = df_ones['atp_inf']),
        dict(range = [float(min(df_ones['atp_sup-atp_inf'])), float(max(df_ones['atp_sup-atp_inf']))],
             label = 'atp_sup-atp_inf', values = df_ones['atp_sup-atp_inf']),
        dict(range=custom_range_AreaD3,
             label='Numer_One_D3', values=df_ones['Numer_One_D3']),
        dict(range=custom_range_AreaD5,
             label='Numer_One_D5', values=df_ones['Numer_One_D5']),
        dict(range=custom_range_AreaD7,
             label='Numer_One_D7', values=df_ones['Numer_One_D7'])
    ])
)
)
#fig.show()

# Save the plot as an HTML file
fig.write_html('parallel_plot_ones.html')

custom_range_AreaD3 = [0, 12]
custom_range_AreaD5 = [0, 12]
custom_range_AreaD7 = [0, 12]

fig = go.Figure(data=
go.Parcoords(
    line_color='blue',
    dimensions = list([
        dict(range = [float(min(df_twos['k_glu'])), float(max(df_twos['k_glu']))],
             constraintrange = [float(min(df_twos['k_glu'])), float(max(df_twos['k_glu']))],
             label = 'k_glu', values = df_area['k_glu']),
        dict(range = [float(min(df_twos['k_aer'])), float(max(df_twos['k_aer']))],
             tickvals = [float(min(df_twos['k_aer'])), float(max(df_twos['k_aer']))],
             label = 'k_aer', values = df_area['k_aer']),
        dict(range = [float(min(df_twos['k_ana'])), float(max(df_twos['k_ana']))],
             tickvals = sorted(df_twos['k_ana'].unique()),
             label = 'k_ana', values = df_twos['k_ana']),
        dict(range = [float(min(df_twos['k_ene'])), float(max(df_twos['k_ene']))],
             label = 'k_ene', values = df_twos['k_ene']),
        dict(range = [float(min(df_twos['k_prolif'])), float(max(df_twos['k_prolif']))],
             label = 'k_prolif', values = df_twos['k_prolif']),
        dict(range = [float(min(df_twos['atp_sup'])), float(max(df_twos['atp_sup']))],
             label = 'atp_sup', values = df_twos['atp_sup']),
        dict(range = [float(min(df_twos['atp_inf'])), float(max(df_twos['atp_inf']))],
             label = 'atp_inf', values = df_twos['atp_inf']),
        dict(range = [float(min(df_twos['atp_sup-atp_inf'])), float(max(df_twos['atp_sup-atp_inf']))],
             label = 'atp_sup-atp_inf', values = df_twos['atp_sup-atp_inf']),
        dict(range=custom_range_AreaD3,
             label='Numer_Two_D0', values=df_twos['Numer_Two_D3']),
        dict(range=custom_range_AreaD5,
             label='Numer_Two_D0', values=df_twos['Numer_Two_D5']),
        dict(range=custom_range_AreaD7,
             label='Numer_Two_D0', values=df_twos['Numer_Two_D7'])
    ])
)
)
#fig.show()

# Save the plot as an HTML file
fig.write_html('parallel_plot_twos.html')

custom_range_AreaD3 = [0, 20]
custom_range_AreaD5 = [0, 20]
custom_range_AreaD7 = [0, 20]

fig = go.Figure(data=
go.Parcoords(
    line_color='blue',
    dimensions = list([
        dict(range = [float(min(df_spheroids['k_glu'])), float(max(df_spheroids['k_glu']))],
             constraintrange = [float(min(df_spheroids['k_glu'])), float(max(df_spheroids['k_glu']))],
             label = 'k_glu', values = df_spheroids['k_glu']),
        dict(range = [float(min(df_spheroids['k_aer'])), float(max(df_spheroids['k_aer']))],
             tickvals = [float(min(df_spheroids['k_aer'])), float(max(df_spheroids['k_aer']))],
             label = 'k_aer', values = df_spheroids['k_aer']),
        dict(range = [float(min(df_spheroids['k_ana'])), float(max(df_spheroids['k_ana']))],
             tickvals = sorted(df_spheroids['k_ana'].unique()),
             label = 'k_ana', values = df_spheroids['k_ana']),
        dict(range = [float(min(df_spheroids['k_ene'])), float(max(df_spheroids['k_ene']))],
             label = 'k_ene', values = df_spheroids['k_ene']),
        dict(range = [float(min(df_spheroids['k_prolif'])), float(max(df_spheroids['k_prolif']))],
             label = 'k_prolif', values = df_spheroids['k_prolif']),
        dict(range = [float(min(df_spheroids['atp_sup'])), float(max(df_spheroids['atp_sup']))],
             label = 'atp_sup', values = df_spheroids['atp_sup']),
        dict(range = [float(min(df_spheroids['atp_inf'])), float(max(df_spheroids['atp_inf']))],
             label = 'atp_inf', values = df_spheroids['atp_inf']),
        dict(range = [float(min(df_spheroids['atp_sup-atp_inf'])), float(max(df_spheroids['atp_sup-atp_inf']))],
             label = 'atp_sup-atp_inf', values = df_spheroids['atp_sup-atp_inf']),
        dict(range=custom_range_AreaD3,
             label='Sp_D3', values=df_spheroids['Sp_D3']),
        dict(range=custom_range_AreaD5,
             label='Sp_D5', values=df_spheroids['Sp_D5']),
        dict(range=custom_range_AreaD7,
             label='Sp_D7', values=df_spheroids['Sp_D7'])
    ])
)
)
#fig.show()

# Save the plot as an HTML file
fig.write_html('parallel_plot_spheroids.html')
