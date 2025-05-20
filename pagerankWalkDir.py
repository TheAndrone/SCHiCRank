import os
import pandas as pd
import networkx as nx
from kneed import KneeLocator
import matplotlib.pyplot as plt
from collections import defaultdict
import pickle

# Configuration
INPUT_DIR = "allResultsWithPhases/K4"
K = 5
MIN_ACTIVE_CELLS = 10

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
    cache_file = os.path.join(input_dir, "full_neighbor_map.pkl")
    if os.path.exists(cache_file):
        print("Loading cached neighbor map...")
        with open(cache_file, "rb") as f:
            return pickle.load(f)
    print("Building full neighbor map...")
    full_neighbor_map = dict()
    for file in os.listdir(input_dir):
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

    with open(cache_file, "wb") as f:
        pickle.dump(full_neighbor_map, f)

    return full_neighbor_map

# Initialize sets
active_cells = set([j for j in range(1056)])
inactive_info = []
iteration = 0

full_neighbor_map = build_full_neighbor_map(INPUT_DIR)

while True:
    cell_pageranks = {cell: [] for cell in active_cells}

    for file, neighbor_dict in full_neighbor_map.items():
        G = nx.DiGraph()
        for cell in active_cells:
            neighbors = [n for n, _ in neighbor_dict.get(cell, []) if n in active_cells][:K]
            for neighbor in neighbors:
                G.add_edge(cell, neighbor)

        # Compute PageRank
        pr = nx.pagerank(G)
        for cell, score in pr.items():
            cell_pageranks[cell].append(score)

    # Process pagerank lists
    pagerank_sums = {}
    for cell, ranks in cell_pageranks.items():
        sorted_ranks = sorted(ranks)
        trimmed = sorted_ranks[2:-2] if len(sorted_ranks) > 4 else sorted_ranks
        pagerank_sums[cell] = sum(trimmed)

    # Sort by sum
    sorted_cells = sorted(pagerank_sums.items(), key=lambda x: x[1], reverse=True)
    values = [score for _, score in sorted_cells]

    # Elbow detection
    x = list(range(len(values)))
    knee_locator = KneeLocator(x, values, curve='convex', direction='decreasing')
    elbow_index = knee_locator.knee

    if elbow_index is None or len(active_cells) <= MIN_ACTIVE_CELLS or elbow_index == len(values):
        break

    # Mark cells to deactivate
    newly_inactive = sorted_cells[elbow_index:]
    for cell, score in newly_inactive:
        inactive_info.append([cell, iteration, score])

    # Update active set
    active_cells = set(cell for cell, _ in sorted_cells[:elbow_index])

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(x, values, marker='o')
    if elbow_index is not None:
        plt.axvline(x=elbow_index, color='red', linestyle='--', label=f'Elbow at {elbow_index}')
    plt.title(f"Iteration {iteration} - PageRank Sum with Elbow")
    plt.xlabel("Cells (sorted)")
    plt.ylabel("Trimmed Sum of PageRanks")
    plt.legend()
    plt.tight_layout()
    plt.show()

    iteration += 1
