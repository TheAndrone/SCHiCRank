import os
import pandas as pd
import networkx as nx
from kneed import KneeLocator
import matplotlib.pyplot as plt
from collections import defaultdict
import pickle


# Configuration
K = 5 # Number of top neighbors to consider for each cell
MIN_ACTIVE_CELLS = 10 # Minimum number of active cells to keep in the analysis

# Load cell phase metadata (global)
metaFn = "cellAndPhaseInfo.pkl"
with open(metaFn, "rb") as f:
    meta = pickle.load(f)
    cell_phases = meta["cell_phase"]


# Helper to get all unique cell IDs from the directory
def get_all_cells(input_dir):
    cell_ids = set()
    for file in os.listdir(input_dir):
        if file.endswith(".csv"):
            df = pd.read_csv(os.path.join(input_dir, file))
            cell_ids.update(df["Item 1"].unique())
            cell_ids.update(df["Item 2"].unique())
    return cell_ids

# Helper to build a master neighbor map for all files
def build_full_neighbor_map(input_dir):
    """
    Builds or loads a full neighbor map from CSV files in the specified input directory.
    This function processes all CSV files in the given directory, where each file contains
    pairwise relationships (e.g., between genomic loci or cells) and their associated frequencies.
    For each item, it constructs a dictionary mapping each entity to its neighbors, sorted by frequency.
    The resulting neighbor maps for all files are stored in a dictionary keyed by filename.
    To improve efficiency, the function caches the result as a pickle file ("full_neighbor_map.pkl")
    in the input directory and loads it on subsequent calls if available.
    Args:
        input_dir (str): Path to the directory containing input CSV files.
    Returns:
        dict: A dictionary where each key is a CSV filename and each value is a defaultdict(list)
              mapping each item to a list of tuples (neighbor, frequency), sorted by descending frequency.
    Notes:
        Each CSV file is expected to have at least the following columns:
        - "Item 1": The first entity in the cell pair (it will be one of the nodes in kNN graph).
        - "Item 2": The second entity in the cell pair.
        - "Frequency": The frequency or strength of the relationship between "Item 1" and "Item 2" 
                        (i.e. number of motiffs two cells share; a metric of similarity)
    """

    cache_file = os.path.join(input_dir, "full_neighbor_map.pkl")
    #If top neighbours for each cells are already computed, just load it
    if os.path.exists(cache_file):
        print("Loading cached neighbor map...")
        with open(cache_file, "rb") as f:
            return pickle.load(f)
        
    #Otherwise, compute the top neighbours for each cell
    print("Building full neighbor map...")
    full_neighbor_map = dict()
    for file in os.listdir(input_dir):
        #File has data for one chromsome
        if not file.endswith(".csv"):
            continue
        print(f"Processing {file}...")
        df = pd.read_csv(os.path.join(input_dir, file))
        neighbor_dict = defaultdict(list)
        for _, row in df.iterrows():
            a, b = row["Item 1"], row["Item 2"]
            freq = row["Frequency"]
            neighbor_dict[a].append((b, freq))
            neighbor_dict[b].append((a, freq))
        for cell in neighbor_dict:
            neighbor_dict[cell].sort(key=lambda x: -x[1])
        full_neighbor_map[file] = neighbor_dict

    #Save for future runs
    with open(cache_file, "wb") as f:
        pickle.dump(full_neighbor_map, f)

    return full_neighbor_map

# Main function

def run_pagerank_filter(INPUT_DIR, label="test", plots=True):
    """
    Iteratively filters cells based on PageRank scores computed from cell k nearest neighbor graphs across chromosomes.
    This function builds directed graphs for each chromosome, where nodes represent cells and edges represent
    top-K neighbors among currently active cells. It computes PageRank scores for each cell, aggregates these
    scores across chromosomes, and removes cells with the lowest aggregate PageRank values in each iteration,
    as determined by an elbow detection algorithm. The process repeats until no further cells can be removed
    or a minimum number of active cells is reached.
    Optionally, the function can plot the PageRank distributions and elbow points for up to three iterations
    at a time.
    Parameters:
        INPUT_DIR (str): Path to the input directory containing data required to build neighbor maps.
        label (str, optional): Label used for output file naming. Default is "test".
        plots (bool, optional): Whether to generate and display plots of PageRank distributions and elbow points.
            Default is True.
    Outputs:
        - Saves a CSV file listing all cells, the iteration in which they were deemed central,
          their final PageRank score, and their phase.
        - Optionally displays plots of PageRank distributions and elbow points for up to three iterations at a time.
    Notes:
        - Requires global variables: `cell_phases`, `K`, `MIN_ACTIVE_CELLS`, and the functions
          `build_full_neighbor_map`, `KneeLocator`, as well as the libraries `networkx`, `matplotlib.pyplot`, and `pandas`.
        - The function assumes that all cells are initially active and iteratively deactivates cells with the lowest
          PageRank scores until the elbow point or a minimum threshold is reached.
    """
    
    resFn = f"final_active_cells_{label}.csv"

    active_cells = set([j for j in range(len(cell_phases))]) #Initially all cells
    inactive_info = []
    iteration = 0

    full_neighbor_map = build_full_neighbor_map(INPUT_DIR)

    while True:
        cell_pageranks = {cell: [] for cell in active_cells}

        for file, neighbor_dict in full_neighbor_map.items():
            #For each chromsome build a directed graph where each cell is a node and edges are the top K neighbours
            G = nx.DiGraph()
            for cell in active_cells:
                neighbors = [n for n, _ in neighbor_dict.get(cell, []) if n in active_cells][:K]
                for neighbor in neighbors:
                    G.add_edge(cell, neighbor)

            # Compute PageRank
            pr = nx.pagerank(G)
            for cell, score in pr.items():
                cell_pageranks[cell].append(score)

        # Process pagerank lists by removing top and bottom 2 values and assigning a mean pagerank to each cell
        pagerank_sums = {}
        for cell, ranks in cell_pageranks.items():
            sorted_ranks = sorted(ranks)
            trimmed = sorted_ranks[2:-2] if len(sorted_ranks) >= 10 else sorted_ranks
            pagerank_sums[cell] = sum(trimmed)

        # Sort by sum
        sorted_cells = sorted(pagerank_sums.items(), key=lambda x: x[1], reverse=True)
        values = [score for _, score in sorted_cells]

        # Elbow detection
        x = list(range(len(values)))
        knee_locator = KneeLocator(x, values, curve='convex', direction='decreasing')
        elbow_index = knee_locator.knee

        # Optionally plot the progress
        if plots:
            # Save iteration data for batch plotting
            if 'batch_plots' not in locals():
                batch_plots = []
            batch_plots.append((iteration+1, x, values, sorted_cells, elbow_index))

            if len(batch_plots) == 3:
                fig, axes = plt.subplots(1, 3, figsize=(18, 5))
                for ax, (it_num, x_vals, y_vals, scells, eidx) in zip(axes, batch_plots):
                    phase_colors = {p: plt.cm.tab10(i % 10) for i, p in enumerate(sorted(set(cell_phases)))}
                    cell_ids = [cell for cell, _ in scells]
                    cell_colors = [phase_colors[cell_phases[cell]] for cell in cell_ids]
                    ax.scatter(x_vals, y_vals, c=cell_colors, s=10)
                    if eidx is not None:
                        ax.axvline(x=eidx, color='red', linestyle='--', label=f'Elkonis pēc {eidx} šūnām')
                    ax.set_title(f"Iterācija {it_num}", fontsize=20)
                    ax.set_xlabel("Šūnas sakārtotas pēc Pagerank", fontsize=18)
                    ax.set_ylabel("PageRank vērtība", fontsize=18)
                handles = [plt.Line2D([0], [0], marker='o', color='w', label=phase,
                                    markerfacecolor=color, markersize=6)
                        for phase, color in phase_colors.items()]
                fig.legend(handles=handles, title="Šūnu fāzes", loc='upper right', fontsize='large', title_fontsize='x-large')
                plt.tight_layout()
                plt.show()
                batch_plots = []

        # If no cell before elbow, finish... could be changed in future to e.g. change k
        if elbow_index is None or len(active_cells) <= MIN_ACTIVE_CELLS or elbow_index == len(values) or elbow_index == 0:
            break

        # Mark cells to deactivate (those on the left of the elbow)
        newly_inactive = sorted_cells[:elbow_index]
        for cell, score in newly_inactive:
            phase = cell_phases[cell] if cell < len(cell_phases) else "Unknown"
            inactive_info.append([cell, iteration, score, phase])

        # Update active set (cells to the right of elbow remain active)
        active_cells = set(cell for cell, _ in sorted_cells[elbow_index:])

        iteration += 1
        print(f"{label} Iteration {iteration}: {len(active_cells)} active cells, elbow at {elbow_index}")

    # Append final active cells to inactive_info
    for cell in sorted(active_cells):
        phase = cell_phases[cell] if cell < len(cell_phases) else "Unknown"
        inactive_info.append([cell, iteration+1, 0, phase])

    # Save inactive cells
    inactive_df = pd.DataFrame(inactive_info, columns=["Cell", "Iteration", "Score", "Phase"])
    inactive_df.to_csv(resFn, index=False)
    print(f"Inactive cells saved to {resFn}")


if __name__ == "__main__":

    run_pagerank_filter("K4_imputed_long_3.0", label="K4_imputed_long_3.0", plots=True) 
