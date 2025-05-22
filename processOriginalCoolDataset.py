from MulticoolProcessor import MulticoolProcessor
from collections import defaultdict
import logging
import os
import json
import pickle


def combineDicts(d1, d2):
    d3 = defaultdict(dict)
    chrs = set(d1.keys()).union(set(d2.keys()))
    for ch in chrs:
        d3[ch] = defaultdict(list)
        
        for (A, B), listOfCells in d1[ch].items():
            if (A, B) in d2[ch]:
                d3[ch][(A, B)] = listOfCells + d2[ch][(A, B)]
            else:
                d3[ch][(A, B)] = listOfCells
        
        for (A, B), listOfCells in d2[ch].items():
            if (A, B) not in d3[ch]:
                d3[ch][(A, B)] = listOfCells
    return d3



def process_cells(cellCount = None, k=1, 
                  chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX"], 
                  fn="sourceData/nagano_10kb_cell_types.scool", 
                  postfix="base10k"):
    
    resolution = 10000
    cellNamesToIndexFn = f"cellNameIndex_{postfix}"
    resFn = f"cellsPerInteraction_{postfix}"
    tmpresFn = f"tmpcellsPerInteraction_{postfix}"

    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    # Initialize the MulticoolProcessor and get cell names
    M = MulticoolProcessor(fn)
    cellNames = M.getCellNames()
    

    # Create a dictionary to map each cell name to an index if not already saved
    cellNamesToIndex = {}
    if os.path.exists(cellNamesToIndexFn):
        with open(cellNamesToIndexFn, 'r') as f:
            cellNamesToIndex = json.load(f)
    else:
        logging.info("Creating cell name index mapping.")
        for idx, cellName in enumerate(cellNames):
            cellNamesToIndex[cellName] = idx
        with open(cellNamesToIndexFn, 'w') as f:
            json.dump(cellNamesToIndex, f)

    if cellCount is not None:
        cellNames = cellNames[:cellCount]

    alreadyProcessedCells = set()
    cellsPerInteractionFull = defaultdict(dict)
    cellsPerInteractionIncrement = defaultdict(dict)

    for i, cellName in enumerate(cellNames):
        if cellName in alreadyProcessedCells:
            logging.info(f"Skipping already processed cell {i}: {cellName}")
            continue  # Skip already processed cells

        logging.info(f"Processing cell {i}: {cellName}")
        try:
            # Process the cell
            cell_index = cellNamesToIndex[cellName]
            C = M.readCell(cellName=cellName)
            if k > 1:
                C = C.reduceResolution(k, fn="tmp.cool")

            for ch in chromosomes:  # Process specific chromosome(s)
                if not C.hasDataOnChr(ch):
                    continue
                chrInteractions = C.getAllInteractionsWithLoci(ch=ch)

                if ch not in cellsPerInteractionIncrement:
                    cellsPerInteractionIncrement[ch] = defaultdict(list)

                for (A, B, count) in chrInteractions:
                    if A==B:
                        continue
                    cellsPerInteractionIncrement[ch][(A, B)].append(cell_index)

            # Mark cell as processed
            alreadyProcessedCells.add(cellName)
        except Exception as e:
            logging.error(f"Error processing cell {cellName}: {e}", exc_info=True)
        
        if ((i+1)%16)==0:
            cellsPerInteractionFull = combineDicts(cellsPerInteractionFull, cellsPerInteractionIncrement)
            cellsPerInteractionIncrement = defaultdict(dict)
            logging.info(f"Processed {i+1} cells.") 
            if ((i+1)%64)==0:
                logging.info("Saving processed cells.")
                with open(f"{tmpresFn}.pkl", 'wb') as f:
                    pickle.dump(cellsPerInteractionFull, f)
                logging.info("Processed cells saved.")

        
        

    cellsPerInteractionFull = combineDicts(cellsPerInteractionFull, cellsPerInteractionIncrement)
    logging.info("Processing completed.")
    logging.info("Saving processed cells.")
    with open(f"{resFn}.pkl", 'wb') as f:
        pickle.dump(cellsPerInteractionFull, f)

    # 
    for cch, dictOfLinkCells in cellsPerInteractionFull.items():
        finalData = {
            "chr": cch,
            "type": postfix,
            "resolution": resolution,
            "cell_links": {},
            "link_cells": {(int(key[0]), int(key[1])): listOfCellsForCurLink for key, listOfCellsForCurLink in dictOfLinkCells.items()},
        }
        #Calculating reverse mapping for cell_links
        for key, listOfCellsForCurLink in dictOfLinkCells.items():
            for cellIndex in listOfCellsForCurLink:
                if cellIndex not in finalData["cell_links"]:
                    finalData["cell_links"][cellIndex] = []
                finalData["cell_links"][cellIndex].append((int(key[0]), int(key[1])))
        #Sort each listOfCellsForCurLink list
        for key, listOfCellsForCurLink in finalData["cell_links"].items():
            finalData["cell_links"][key] = sorted(listOfCellsForCurLink)
        
        #Save to pkl
        finalFn = f"{postfix}-{cch}-{resolution}.pkl"

        with open(finalFn, 'wb') as f:
            pickle.dump(finalData, f)


if __name__ == "__main__":
    process_cells(cellCount=None, k=1, 
                  chromosomes=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX"],
                  fn="sourceData/nagano_10kb_cell_types.scool", 
                  postfix="base10k")

