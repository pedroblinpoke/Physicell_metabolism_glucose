import pandas as pd
import matplotlib.pyplot as plt

# Load the data from a .txt file (assuming space or tab-separated)
file_path = 'mean_meta_area_no12.txt'  # Replace with your file's path
df = pd.read_csv(file_path, delim_whitespace=True)

# Check if the file is loaded correctly
print(df.head())

# Plotting Area7 vs Metabolism with error bars using Var7
plt.errorbar(df['Metabolism'], df['Mean_AreaD7no12'], yerr=df['SEM_AreaD7no12'], fmt='o', capsize=5, label='Area7 with Var7 error')


plt.xscale('log')
# Labeling the plot
plt.xlabel('Metabolism')
plt.ylabel('Area')
plt.title('Area vs Metabolism Dia 7 2D')
plt.legend()
plt.grid(True)
plt.show()
