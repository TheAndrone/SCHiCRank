# SCHiCRank
## Overview

**SCHiCRank** is an iterative cell filtering algorithm that leverages PageRank centrality on k-nearest neighbor (kNN) graphs, constructed from pairwise cell similarity data across chromosomes.

---

## Dependencies

- `pandas`
- `networkx`
- `kneed`
- `matplotlib`
- `collections` (`defaultdict`)
- `pickle`
- `os`

---

## Workflow

1. **Load cell phase metadata** from a pickle file.
2. **Build or load a neighbor map** from CSV files in the input directory, where each file contains pairwise cell relationships and similarity scores.
3. **Iteratively construct kNN graphs** for each chromosome using active cells, compute PageRank scores, and aggregate these scores across chromosomes.
4. **Identify the elbow point** in sorted aggregate PageRank scores using the KneeLocator algorithm; deactivate cells with the lowest scores (left of the elbow).
5. **Repeat** until no further cells can be removed or a minimum number of active cells is reached.
6. **Optionally plot** PageRank distributions and elbow points for up to three iterations at a time.
7. **Save results** as a CSV file listing all cells, their iteration of centrality, final PageRank score, and phase.

---

## Functions

- `get_all_cells(input_dir)`: Returns all unique cell IDs from the input directory.
- `build_full_neighbor_map(input_dir)`: Builds or loads a neighbor map for all files in the input directory, mapping each cell to its sorted neighbors by frequency.
- `run_pagerank_filter(INPUT_DIR, label="test", plots=True)`: Main function that performs iterative PageRank-based filtering and outputs results.

---

## Usage

Run the script directly to execute the PageRank filtering on the specified input directory.

