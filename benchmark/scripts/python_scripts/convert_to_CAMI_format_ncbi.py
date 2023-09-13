"""
This script is used to converet the YACHT results (based on GTDB reference database) to CAMI profiling Bioboxes format(https://github.com/bioboxes/rfc/blob/60263f34c57bc4137deeceec4c68a7f9f810f6a5/data-format/profiling.mkd)
"""

import os
from pathlib import Path
import argparse
import pandas as pd
import requests
import tarfile
import pytaxonkit
import numpy as np
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is used to converet the YACHT results (based on GTDB reference database) to CAMI profiling Bioboxes format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--yacht_res", type=str, required=True, help="path to the YACHT CSV result")
    parser.add_argument("--metadata_dir", type=str, required=False, help="path to the metadata directory", default=None)
    parser.add_argument('--min_coverage', type=float, help='To compute false negative weight, assume each organism has this minimum coverage in sample. Should be between 0 and 1.', required=False, default = 1)
    parser.add_argument("--outfile", type=str, required=True, help="path to the CAMI formatted output file")
    args = parser.parse_args()

    ## allowable_rank
    allowable_rank = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

    ## Read GTDB metadata if available
    if args.metadata_dir is not None:
        ncbi_metadata_path = Path(args.metadata_dir)
        if not ncbi_metadata_path.is_dir():
            print(f"Error: {args.metadata_dir} is not a directory", flush=True)
            exit(1)
        else:
            # Read the metadata
            ncbi_metadata_path = str(ncbi_metadata_path.absolute())
            metadata_df = pd.read_csv(os.path.join(ncbi_metadata_path, 'ncbi_database_metadata.tsv'), sep='\t', header=0)
    else:
        print("Error: NCBI metadata is not provided")
        exit(1)

    ## Find the ncbi lineage
    result = pytaxonkit.lineage(list(set(metadata_df['NCBI_ID'].tolist())))
    metadata_df = metadata_df.merge(result[['TaxID','Rank','FullLineageTaxIDs','FullLineage','FullLineageRanks']], left_on='NCBI_ID', right_on='TaxID').drop(columns=['NCBI_ID']).reset_index(drop=True)
    metadata_df['FullLineageTaxIDs'] = metadata_df['FullLineageTaxIDs'].str.replace(';','|')
    metadata_df['FullLineage'] = metadata_df['FullLineage'].str.replace(';','|')
    metadata_df['FullLineageRanks'] = metadata_df['FullLineageRanks'].str.replace(';','|')

    # Read the YACHT results
    yacht_res_path = Path(args.yacht_res)
    if not yacht_res_path.is_file():
        print(f"Error: {yacht_res_path} is not a file")
        exit(1)
    else:
        yacht_res_path = str(yacht_res_path.absolute())
    yacht_res_df = pd.read_csv(yacht_res_path, sep=',', header=0)
    ## drop the first column
    yacht_res_df.drop(yacht_res_df.columns[0], axis=1, inplace=True)
    ## re-determine the in_sample_est based on the given min_coverage
    yacht_res_df['min_coverage'] = args.min_coverage
    yacht_res_df['in_sample_est'] = yacht_res_df.apply(lambda row: (row['num_matches'] >= row['acceptance_threshold_wo_coverage'] * row['min_coverage']) & (row['num_matches'] != 0) & (row['acceptance_threshold_wo_coverage'] != 0), axis=1)
    ## select the organisms that YACHT considers to present in the sample
    organism_name_list = yacht_res_df.query('in_sample_est == True')['organism_name'].tolist()
    organism_gtdbid_list = [x.split(' ')[0] for x in organism_name_list]

    if len(organism_gtdbid_list) == 0:
        print(f"Error: No organism is detected by YACHT")
        exit(1)

    ## Merge the YACHT results with the metadata
    selected_organism_metadata_df = metadata_df.query("genome_ID in @organism_gtdbid_list").reset_index(drop=True)

    ## Summarize the results
    summary_dict = {}
    for row in selected_organism_metadata_df.to_numpy():
        select_index = [index for index, x in enumerate(row[-1].split('|')) if x in allowable_rank]
        taxid_list = list(np.array(row[5].split('|'))[select_index])
        lineage_list = list(np.array(row[6].split('|'))[select_index])
        rank_list = list(np.array(row[7].split('|'))[select_index])
        current_lineage = ''
        current_taxid = ''
        for index, (taxid, rank, lineage) in enumerate(zip(taxid_list, rank_list, lineage_list)):
            if index == 0:
                current_lineage = lineage
                current_taxid = taxid
            else: 
                current_lineage = current_lineage + '|' + lineage
                current_taxid = current_taxid + '|' + taxid
            if taxid not in summary_dict:
                summary_dict[taxid] = {'RANK':rank, 'TAXPATH':current_taxid, 'TAXPATHSN':current_lineage, 'count':1}
            else:
                summary_dict[taxid]['count'] += 1

    # calculate percentage
    for taxid in summary_dict:
        summary_dict[taxid]['PERCENTAGE'] = summary_dict[taxid]['count']/len(selected_organism_metadata_df)*100

    ## sort by rank in allowable rank list
    summary_df = pd.DataFrame(summary_dict).T.reset_index().rename(columns={'index':'TAXID'})
    res_df = [summary_df.query(f'RANK == "{rank}"') for rank in allowable_rank]
    res_df = pd.concat(res_df).drop(columns=['count']).reset_index(drop=True)
    res_df.columns = ['@@TAXID', 'RANK', 'TAXPATH', 'TAXPATHSN', 'PERCENTAGE']

    ## output summary results
    sample_name = os.path.basename(yacht_res_path).split('.')[0]
    with open(args.outfile, 'w') as f:
        f.write(f'@SampleID:{sample_name}\n')
        f.write('@Version:0.9.1\n')
        f.write(f"@Ranks:{'|'.join(list(res_df['RANK'].unique()))}\n\n")
        f.write(f"@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
        for row in res_df.to_numpy():
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\n")