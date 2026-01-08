"""
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
from pathlib import Path
import numpy as np
import pandas as pd




file_path_area =  Path(f'D:/Resultados_simulaciones/30_5_2024_Col_25/output_areas_modified.txt')
df_area = pd.read_csv(file_path_area, delim_whitespace=True)

problem = {
    'num_vars': 7,
    'names': ['k_glu', 'k_aer', 'k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf'],
    'bounds': [[0.0000001, 0.001],
               [0.0000001, 0.001],
               [0.0000001, 0.001],
               [0.0000001, 0.001],
               [0.00045, 0.00077],
               [10, 100000],
               [1, 10000]]
}

param_values = df_area.iloc[:, :7]
Y = df_area.iloc[:, -1:]


param_values = param_values.to_numpy()


Y = Y.to_numpy()
Y = Y.flatten()


problem = {
    'num_vars': 3,
    'names': ['x1', 'x2', 'x3'],
    'bounds': [[-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359]]
}

param_values = saltelli.sample(problem, 1024)


Y = Ishigami.evaluate(param_values)
print(f"Params values{param_values}")
print(f"Results{Y}")
print(param_values.shape)    
print(Y.size)


Si = sobol.analyze(problem, Y,calc_second_order=False)
print(f"Sobol indexes:{Si['S1']}")

"""



import pandas as pd

# Load the data from the text file
file_path = 'D:/Resultados_simulaciones/30_5_2024_Col_25/output_areas_modified.txt'  # Replace with the actual path to your file
data = pd.read_csv(file_path, sep='\t')  # Adjust the separator if necessary

# Separate parameters and output
parameters = data[['k_glu', 'k_aer', 'k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf']].values
output = data['AreaD7'].values


from SALib.analyze import sobol

# Define the model input parameters
problem = {
    'num_vars': 7,
    'names': ['k_glu', 'k_aer', 'k_ana', 'k_ene', 'k_prolif', 'atp_sup', 'atp_inf'],
    'bounds': [[data['k_glu'].min(), data['k_glu'].max()],
               [data['k_aer'].min(), data['k_aer'].max()],
               [data['k_ana'].min(), data['k_ana'].max()],
               [data['k_ene'].min(), data['k_ene'].max()],
               [data['k_prolif'].min(), data['k_prolif'].max()],
               [data['atp_sup'].min(), data['atp_sup'].max()],
               [data['atp_inf'].min(), data['atp_inf'].max()]]
}


# Perform the Sobol sensitivity analysis
sobol_indices = sobol.analyze(problem, output, print_to_console=True)

# Print the results
print('Sobol First Order Indices:', sobol_indices['S1'])
print('Sobol Total Order Indices:', sobol_indices['ST'])
print('Sobol Second Order Indices:', sobol_indices['S2'])
