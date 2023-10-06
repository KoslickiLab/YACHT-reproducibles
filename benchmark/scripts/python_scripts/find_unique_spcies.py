"""
This script is used to find the unique species-level genome. That is, if there are multiple strain genomes under the same species, randomly pick one
"""


import os
from pathlib import Path
import argparse
import pandas as pd
import pytaxonkit
from tqdm import tqdm

def custom_func(group):
    if group['target'].any():  # If any target==True in the group
        return group[group['target']]  # Return all rows with target==True
    else:
        return group.head(1)  # Return only the first row

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is used to find the unique species-level genome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--taxids", type=str, required=True, help="path to a file containing all taxids that will be used for database build")
    parser.add_argument("--target_taxids", type=str, required=False, help="path to a file containing target taxids to keep", default=None)
    parser.add_argument("--output", type=str, required=False, help="path to the output TSV file", default="unique_species_taxids.txt")
    args = parser.parse_args()
    
    ## read in taxids
    taxids = pd.read_csv(args.taxids, header=None, names=["taxid"])
    
    ## find the lineage of each taxid
    lineage_result = pytaxonkit.lineage(taxids['taxid'].to_list())
    
    ## find the species-level taxid of each given taxid
    for index, row in enumerate(tqdm(lineage_result.to_numpy())):
        try:
            species_level_index = row[-1].split(';').index('species')
        except ValueError:
            species_level_index = -1
        if species_level_index != -1:
            species_taxid = row[0]
            lineage_result.loc[index, 'species_taxid'] = row[-2].split(';')[species_level_index]
    
    ## add target taxids if provided
    if args.target_taxids is not None:
        target_taxids = pd.read_csv(args.target_taxids, header=None, names=["taxid"])
        lineage_result['target'] = False
        lineage_result.loc[lineage_result['TaxID'].isin(target_taxids['taxid']), 'target'] = True
    else:
        lineage_result['target'] = False
    
    # find the unique species-level taxid
    filtered_lineage_result = lineage_result.groupby('species_taxid').apply(custom_func).reset_index(drop=True)
    
    # write to file
    filtered_lineage_result.to_csv(args.output, sep='\t', index=None)