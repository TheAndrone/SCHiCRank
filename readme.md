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

Download the `nagano_10kb_cell_types.scool` and `nagano_assoziated_cell_types.txt` files from [Zenodo](https://zenodo.org/records/4308298) and save them to the `sourceData` directory.

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
    "index_to_name": {<cell_index>: <cell_name>},
    "index_to_type": {<cell_index>: <cell_type_in_mitotic_cycle>},
    "cell_IDs": [<cell_index>],
    "cell_links": {<cell_index>: [(A, B), ...], ...},
    "link_cells": {(A, B): [<cell_index>, ...], ...},
}
```

### Parameters

- **fn** (`str`, optional):  
  Path to the `.scool` file. Default: `"sourceData/nagano_10kb_cell_types.scool"`.

- **fnResolution** (`int`, optional):  
  Resolution of the original `.scool` file. Default: `10000`.
  
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



# Script: [`createCliqueDatafiles.py`](./createCliqueDatafiles.py)


This script processes the files created by `process_cells` and generates new `.pkl` files containing clique information for each cell. 

- **cell_cliques**: For each cell and clique size (K3–K8), a set of cliques (as sorted tuples of node indices) found in that cell.
- **clique_cells**: For each clique size (K3–K8), a mapping from each clique (tuple) to the set of cells where it appears.
- Other metadata fields:  
  - `chr`: Chromosome name  
  - `resolution`: Bin resolution  
  - `type`: Data type (with `_cliques` appended)  
  - `index_to_name`: Mapping from node index to name  
  - `index_to_type`: Mapping from node index to type  
  - `cell_IDs`: List of cell IDs

**Output `.pkl` structure:**
```python
{
    "chr": ...,
    "resolution": ...,
    "type": ...,
    "index_to_name": ...,
    "index_to_type": ...,
    "cell_IDs": [...],
    "cell_cliques": {
        "K3": {cellID: set(cliques), ...},
        ...,
        "K8": {cellID: set(cliques), ...}
    },
    "clique_cells": {
        "K3": {clique: set(cellIDs), ...},
        ...,
        "K8": {clique: set(cellIDs), ...}
    }
}
```
This script adds cliques for each cell and the reverse mapping of cells for each clique. Can be used for further analysis.



# Script: [`createCliqueCountsOverview.py`](./createCliqueCountsOverview.py)


## Overview

This script processes clique data stored in a pickle file and generates a CSV summary of clique counts for each cell. It reads a `.pkl` file containing information about cliques detected in single-cell Hi-C data, extracts relevant statistics (such as the number of cliques of different sizes and lengths per cell), and writes the results to a CSV file.

### Input

- **Pickle file**: Contains a dictionary with clique information, cell metadata, and chromosome details, as generated by [`createCliqueDatafiles.py`](./createCliqueDatafiles.py).

### Output

- **CSV file**: Each row corresponds to a cell, with columns for:
    - `cellID`: Unique identifier for the cell.
    - `cellName`: Name of the cell.
    - `cellPhase`: Cell phase in the mitotic cycle, i.e. the true cluster.
    - `chr`: Chromosome.
    - Additional columns for each clique type and size (e.g., `base100k_cliques-K3-alllengths`, `base100k_cliques-K3-long`, etc.), representing the count of cliques per cell.

The resulting CSV can be used for cell cycle prediction where each cell is characterized by the number of cliques it contains. 




# Script: [`createPairwiseSimilarities.py`](./createPairwiseSimilarities.py) 

This script processes clique data from a pickle file generated by [`createCliqueDatafiles.py`](./createCliqueDatafiles.py) and computes pairwise similarity between each pair of cells based on co-occurrence of cliques in these cells. It outputs CSV files containing the frequency of each cell pair, along with metadata such as cell names, phases, and types. The results are organized by motif and clique length, and are saved in dedicated directories for further analysis.

## Output Format

Each output CSV file contains the following columns:

- **Item 1**: Index of the first cell in the pair.
- **Item 2**: Index of the second cell in the pair.
- **Frequency**: Number of cliques in which this cell pair co-occurs.
- **cell1_name**: Name of the first cell.
- **cell2_name**: Name of the second cell.
- **cell1_phase**: Phase/type annotation for the first cell.
- **cell2_phase**: Phase/type annotation for the second cell.
- **same?**: Boolean indicating if both cells have the same phase/type.
- **type**: Concatenation of both cells' types (e.g., "G1+early-S").

Each CSV is saved in a subdirectory named after the motif and clique length, with filenames encoding the analysis parameters.
"""