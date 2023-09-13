"""
This script is used to build ncbi database metadata file.
"""


import os
from pathlib import Path
import argparse
import pandas as pd
import pytaxonkit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is used to build ncbi database metadata file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", nargs='+', type=str, required=True, help="paths to the 'combined_metadata.tsv' files")
    parser.add_argument("--refseq_taxid_mapping", type=str, required=True, help="path to the 'refseq_taxid_mapping.csv' file")
    args = parser.parse_args()
    
    ## read metadata files
    if len(args.metadata) > 0:
        metadata_path_list = [Path(path) for path in args.metadata]
        for path in metadata_path_list:
            if not path.exists():
                raise FileNotFoundError(f"File '{str(path)}' does not exist.")
        df_list = [pd.read_csv(str(path), sep='\t', header=0).drop(columns=['OTU']) for path in metadata_path_list]
        pos_genome_df = pd.concat(df_list).reset_index(drop=True)
        
    ## read refseq_taxid_mapping file
    refseq_taxid_mapping_path = Path(args.refseq_taxid_mapping)
    if not refseq_taxid_mapping_path.exists():
        raise FileNotFoundError(f"File '{str(refseq_taxid_mapping_path)}' does not exist.")
    refseq_taxid_mapping_df = pd.read_csv(str(refseq_taxid_mapping_path), sep=',', header=None)
    refseq_taxid_mapping_df.columns = ['refseq', 'taxid']
    # remove the NULL row
    refseq_taxid_mapping_df = refseq_taxid_mapping_df.loc[~refseq_taxid_mapping_df['taxid'].isna(),:].reset_index(drop=True)
    refseq_taxid_mapping_df['taxid'] = refseq_taxid_mapping_df['taxid'].astype(int)
    
    ## Find the ncbi lineage
    result = pytaxonkit.lineage(list(set(refseq_taxid_mapping_df['taxid'].tolist())))
    chosen_taxids = result.loc[result['Rank'].isin(['species', 'strain']),'TaxID'].to_list()
    refseq_taxid_mapping_df = refseq_taxid_mapping_df.loc[refseq_taxid_mapping_df['taxid'].isin(chosen_taxids),:].reset_index(drop=True)
    ## remove the species that are in the pos_genome_df
    refseq_taxid_mapping_df = refseq_taxid_mapping_df.loc[~refseq_taxid_mapping_df['taxid'].isin(pos_genome_df['NCBI_ID'].tolist()),:].reset_index(drop=True)
    refseq_taxid_mapping_df['novelty_category'] = 'known_ncbi'
    temp_df = pd.DataFrame([('_'.join(os.path.basename(x).split('_')[:2]), str(x)) for x in Path(os.path.join(os.path.dirname(refseq_taxid_mapping_path.absolute()), 'RefSeq_genomic_20190108')).glob('*.fna.gz')])
    temp_df.columns = ['refseq', 'path']
    refseq_taxid_mapping_df = refseq_taxid_mapping_df.merge(temp_df, on='refseq').reset_index(drop=True)
    refseq_taxid_mapping_df.columns = ['genome_ID', 'NCBI_ID', 'novelty_category', 'path']
    
    ## Combine two tables
    final_df = pd.concat([pos_genome_df, refseq_taxid_mapping_df]).reset_index(drop=True)
    final_df['genome_ID'] = [f'genome{i+1}'for i in range(final_df.shape[0])]
    final_df.to_csv(str(Path(os.path.dirname(refseq_taxid_mapping_path.absolute()), 'ncbi_database_metadata.tsv')), sep='\t', index=None)

