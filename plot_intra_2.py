from typing import Dict, List
from scipy import io as sio
from pathlib import Path
import pandas as pd
import physicool.processing as processing
import seaborn as sns
import matplotlib.pyplot as plt
from xml.etree import ElementTree

# Define where to find the output files
OUTPUT_PATH = Path("output")
# Define the cell variables to be extracted (the script has a method to search for their indices)
# VARIABLES = ["intra_oxy", "intra_glu", "intra_lac", "intra_energy", "intra_pyr", "intra_NAD", "intra_NADH"]
VARIABLES = ["intra_oxy", "intra_glu", "intra_lac", "intra_energy", "intra_NAD", "intra_NADH", "intra_pyr", "k_gly", "k_e", "total_volume"]


def get_variables_idx(variables: List[str], output_path: Path) -> Dict[str, int]:
    """Gets the custom data variables and their indices from the initial.xml file."""
    tree = ElementTree.parse(output_path / "initial.xml")
    stem = "cellular_information/cell_populations/cell_population/custom/simplified_data[@source='PhysiCell']/labels"
    return {label.text: int(label.attrib["index"]) for label in tree.find(stem).findall("label") if
            label.text in variables}


def extract_data(frame: int, variables_map: Dict[str, int], output_path: Path) -> pd.DataFrame:
    """Extracts the passed variables from the output files for the given time point."""
    time_str = str(frame).zfill(8)
    path_name = output_path / f"output{time_str}_cells_physicell.mat"
    cells = sio.loadmat(path_name)["cells"]
    # Create a DataFrame where each column is a cell variable
    data = pd.DataFrame({key: cells[value] for key, value in variables_map.items()})
    # Add the current time point to the DataFrame to identify it later
    data["time"] = frame
    return data


if __name__ == "__main__":
    # Get the indices of the variables to be extracted
    variables_map = get_variables_idx(output_path=OUTPUT_PATH, variables=VARIABLES)
    print(variables_map)
    # Get the number of output files in the output folder
    frames = processing.get_cell_file_num(output_path=OUTPUT_PATH, version="1.9.0")
    # Create a DataFrame with all the data for all the time points
    df = pd.concat([extract_data(frame, variables_map, output_path=OUTPUT_PATH) for frame in range(frames)],
                   ignore_index=True)

    # df.to_excel("substances_intra.xlsx")

    # Plot the results
    fig, axes = plt.subplots(2, 5, figsize=(14, 4))
    # The ravel() function "flattens" the axes array from two columns
    # (e.g., 2x2) to one column (e.g., 1x4) so that only one loop is needed.
    # enumerate() will give the ax object and its index in the flattened array.
    # We can use this index to access the variables in the VARIABLES list
    # (each plot will be a different variable, according to his index).
    for i, ax in enumerate(axes.ravel()):
        sns.lineplot(data=df, x="time", y=VARIABLES[i], color=f"C{i}", ax=ax)
        ax.set_title(f"Substance: {VARIABLES[i]}")

    sns.despine()
    plt.tight_layout()
    plt.show()
