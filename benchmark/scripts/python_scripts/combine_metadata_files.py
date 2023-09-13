# sourcery skip: collection-builtin-to-comprehension, comprehension-to-generator, for-index-underscore
"""
This script is used to combine the metadata files.
"""


import os
from pathlib import Path
import argparse
import pandas as pd
import pytaxonkit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is used to combine the metadata files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", type=str, required=True, help="path to 'metadata.tsv' file")
    parser.add_argument("--genome_to_id", type=str, required=True, help="path to the 'genome_to_id.tsv' file")
    parser.add_argument("--ground_truth_dir", type=str, required=True, help="path to the ground truth directory")
    args = parser.parse_args()
    
    ## Read the metadata file
    metadata_path = Path(args.metadata)
    if not metadata_path.is_file():
        raise Exception(f"ERROR: The metadata file '{args.metadata}' does not exist.")
    metadata_df = pd.read_csv(metadata_path, sep="\t", header=0)
    
    ## Read the genome_to_id file
    genome_to_id_path = Path(args.genome_to_id)
    if not genome_to_id_path.is_file():
        raise Exception(f"ERROR: The genome_to_id file '{args.genome_to_id}' does not exist.")
    genome_to_id_df = pd.read_csv(genome_to_id_path, sep="\t", header=None)
    genome_to_id_df[1] = genome_to_id_df[1].apply(lambda cell: os.path.join(os.path.dirname(metadata_path.absolute()), os.path.basename(cell)))
    genome_to_id_df.columns = ["genome_ID", "path"]
    
    ## merge two tables
    combined_df = pd.merge(metadata_df, genome_to_id_df, on="genome_ID").reset_index(drop=True)

    ## Read the ground truth files
    ground_truth_dir = Path(args.ground_truth_dir)
    if not ground_truth_dir.is_dir():
        raise Exception(f"ERROR: The ground truth directory '{args.ground_truth_dir}' does not exist.")
    ground_truth_df = [pd.read_csv(x, sep='\t', skiprows=3) for x in ground_truth_dir.glob("*_sample_*.profile") if 'long' not in str(x)]
    ## Get all taxoids from the ground truth files
    taxid_list = list(set([int(taxid) for df in ground_truth_df for taxid in df['@@TAXID'].to_list()]))
    combined_df = combined_df.query('NCBI_ID in @taxid_list').reset_index(drop=True)
    ## remove unidentified species
    result = pytaxonkit.lineage(set(combined_df['NCBI_ID']))
    combined_df = combined_df.loc[~combined_df['NCBI_ID'].isin(result.loc[result['Name'].str.contains('unidentified'),'TaxID'].to_list()),:]
    combined_df.to_csv(os.path.join(os.path.dirname(metadata_path.absolute()), "combined_metadata.tsv"), sep="\t", index=False)