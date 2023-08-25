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

def download_file(url, save_path):
    response = requests.get(url, stream=True)

    # Check if the request was successful
    if response.status_code == 200:
        # Write the response content to a local file
        with open(os.path.join(save_path), 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    else:
        print(f"Failed to download the file. HTTP status code: {response.status_code}")

def read_tar_gz_file(tar_gz_file_path):
    tar = tarfile.open(tar_gz_file_path, "r:gz")
    return pd.read_csv(tar.extractfile(tar.getmembers()[0]), sep='\t', header=0)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is used to converet the YACHT results (based on GTDB reference database) to CAMI profiling Bioboxes format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--yacht_res", type=str, required=True, help="path to the YACHT CSV result")
    parser.add_argument("--gtdb_metadata_dir", type=str, required=False, help="path to the GTDB metadata directory", default=None)
    parser.add_argument("--outfile", type=str, required=True, help="path to the CAMI formatted output file")
    args = parser.parse_args()

    ## allowable_rank
    allowable_rank = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

    ## Read GTDB metadata if available
    if args.gtdb_metadata_dir is not None:
        gtdb_metadata_path = Path(args.gtdb_metadata_dir)
        if not gtdb_metadata_path.is_dir():
            print(f"Error: {args.gtdb_metadata_dir} is not a directory", flush=True)
            exit(1)
        elif len(list(gtdb_metadata_path.glob('*.tar.gz'))) != 2:
            print(f"Error: {args.gtdb_metadata_dir} does not contain any '*.tar.gz' file. Download it automaticallyy.", flush=True)
            download_file("https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_metadata_r214.tar.gz", os.path.join(args.gtdb_metadata_dir,'ar53_metadata_r214.tar.gz'))
            download_file("https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_metadata_r214.tar.gz", os.path.join(args.gtdb_metadata_dir,'bac120_metadata_r214.tar.gz'))
            # Read the metadata
            gtdb_metadata_path = str(gtdb_metadata_path.absolute())
            bac120_metadata_df = read_tar_gz_file(os.path.join(gtdb_metadata_path,'bac120_metadata_r214.tar.gz'))[['accession','ncbi_taxid']]
            ar53_metadata_df = read_tar_gz_file(os.path.join(gtdb_metadata_path,'ar53_metadata_r214.tar.gz'))[['accession','ncbi_taxid']]
            metadata_df = pd.concat([ar53_metadata_df, bac120_metadata_df]).reset_index(drop=True)
            metadata_df['accession'] = metadata_df['accession'].str.replace('^(GB_|RS_)', '', regex=True)
        else:
            # Read the metadata
            gtdb_metadata_path = str(gtdb_metadata_path.absolute())
            bac120_metadata_df = read_tar_gz_file(os.path.join(gtdb_metadata_path,'bac120_metadata_r214.tar.gz'))[['accession','ncbi_taxid']]
            ar53_metadata_df = read_tar_gz_file(os.path.join(gtdb_metadata_path,'ar53_metadata_r214.tar.gz'))[['accession','ncbi_taxid']]
            metadata_df = pd.concat([ar53_metadata_df, bac120_metadata_df]).reset_index(drop=True)
            metadata_df['accession'] = metadata_df['accession'].str.replace('^(GB_|RS_)', '', regex=True)
    else:
        print("Error: GTDB metadata is not provided")
        exit(1)

    ## Find the ncbi lineage
    result = pytaxonkit.lineage(list(set(metadata_df['ncbi_taxid'].tolist())))
    metadata_df = metadata_df.merge(result[['TaxID','Rank','FullLineageTaxIDs','FullLineage','FullLineageRanks']], left_on='ncbi_taxid', right_on='TaxID').drop(columns=['ncbi_taxid']).reset_index(drop=True)
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
    ## select the organisms that YACHT considers to present in the sample
    organism_name_list = yacht_res_df.query('in_sample_est == True')['organism_name'].tolist()
    organism_gtdbid_list = [x.split(' ')[0] for x in organism_name_list]

    if len(organism_gtdbid_list) == 0:
        print(f"Error: No organism is detected by YACHT")
        exit(1)

    ## Merge the YACHT results with the metadata
    selected_organism_metadata_df = metadata_df.query("accession in @organism_gtdbid_list").reset_index(drop=True)

    ## Summarize the results
    summary_dict = {}
    for row in selected_organism_metadata_df.to_numpy():
        select_index = [index for index, x in enumerate(row[-1].split('|')) if x in allowable_rank]
        taxid_list = list(np.array(row[3].split('|'))[select_index])
        lineage_list = list(np.array(row[4].split('|'))[select_index])
        rank_list = list(np.array(row[5].split('|'))[select_index])
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
        f.write(f"{'|'.join(list(res_df['RANK'].unique()))}\n")
        for row in res_df.to_numpy():
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\n")