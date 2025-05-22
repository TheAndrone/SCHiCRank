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

Run the script [`runSCHiCRank.py`](./runSCHiCRank.py) directly to execute the SCHiCRank on an example directory.




# Pipeline to process original data in scool format for clique-based topological analysis
## Function: `process_cells` in [`processOriginalCoolDataset.py`](./processOriginalCoolDataset.py)

The `process_cells` function processes a `.scool` file containing scHi-C data for multiple cells, extracting and saving interactions for specified chromosomes into a format suitable for downstream topological analysis.

### Download Example Data

Download the `nagano_10kb_cell_types.scool` file from [Zenodo](https://zenodo.org/records/4308298) and save it to the `sourceData` directory.

# Reference
This file is from "Robust and efficient single-cell Hi-C clustering with approximate k-nearest neighbor graphs".

**Creators:** Joachim Wolff

**Description:**  
Date used in 'Robust and efficient single-cell Hi-C clustering with approximate k-nearest neighbor graphs'.

**Original data source:**  
- Nagano 2017: GEO94489  


### Output Format

```python
{
    "chr": <chromosome name>,
    "type": <identifier, e.g., "base10k">,
    "resolution": <resolution in base pairs, e.g., 10000>,
    "cell_links": {<cell_index>: [(A, B), ...], ...},
    "link_cells": {(A, B): [<cell_index>, ...], ...},
}
```

### Parameters

- **fn** (`str`, optional):  
  Path to the `.scool` file. Default: `"sourceData/nagano_10kb_cell_types.scool"`.

- **postfix** (`str`, optional):  
  Postfix for output files. Default: `"base10k"`.

- **chromosomes** (`list of str`, optional):  
  List of chromosomes to process. Default: `["chr1", ..., "chrX"]`.

- **cellCount** (`int`, optional):  
  Number of cells to process. If `None`, all cells are processed. Default: `None`.

- **k** (`int`, optional):  
  Coarsening factor for reducing resolution. Default: `1` (no reduction).
    
---

### Additional Dependencies

- `h5py`
- `cooler`

Helper functions used in `process_cells` are implemented in the `MulticoolProcessor` and `CoolProcessor` classes.
