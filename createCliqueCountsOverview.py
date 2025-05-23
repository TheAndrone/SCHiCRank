import pickle
import os


# Helper function to process cliques
def process_cliques(fn, resFN):

    counts = {}
    print(f"Processing cliques: {fn}")

    with open(fn, "rb") as f:
        data = pickle.load(f)

    typeName = data["type"]
    for K in data["cell_cliques"].keys():
        #K is like K3, K4,...
        for cliqueSize in ["alllengths", "long"]:
            finalLabel = f"{typeName}-{K}-{cliqueSize}"
            print(f"Processing cliques for {finalLabel}...")
            cliques = data["cell_cliques"][K]
            for cellID, cliques in cliques.items():
                if cellID not in counts:
                    counts[cellID] = {"cellID": cellID, "cellName": data["index_to_name"][cellID], "cellPhase": data["index_to_type"][cellID], "chr": data["chr"]}
                eligibleCliques = len([clique for clique in cliques if clique[-1]-clique[0]>=2000000]) if cliqueSize == "long" else len(cliques)
                counts[cellID][finalLabel] = eligibleCliques 
            
    save_counts(counts, resFN)


def save_counts(countsD, resFN):
    if not countsD:
        print("No data found")
        return
    
    # Collect all possible keys across all cells
    all_keys = set()
    for cell_data in countsD.values():
        all_keys.update(cell_data.keys())

    # Prioritize standard keys
    prefix_keys = ["cellID", "cellName", "cellPhase", "chr"]
    other_keys = sorted(k for k in all_keys if k not in prefix_keys)
    ordered_keys = prefix_keys + other_keys

    # Write to CSV
    with open(resFN, "w") as f:
        f.write(",".join(ordered_keys) + "\n")
        for data in countsD.values():
            row = [str(data.get(k, "")) for k in ordered_keys]
            f.write(",".join(row) + "\n")

    print(f"Saved counts to {resFN}")



if __name__ == "__main__":
    process_cliques("base100k-chr18-100000-cliques.pkl", "all_cliques_data_chr18_base100k.csv")
