# sourcery skip: raise-specific-error
"""
This script is used to generate the plots for comparing YACHT with other SOTA tools.
"""


import os
from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is used to generate the plots for comparing YACHT with other SOTA tools.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cami_opal_res", type=str, required=True, help="path to the CAMI reported OPAL result file")
    parser.add_argument("--yacht_opal_dir", type=str, required=True, help="path to the YACHT OPAL result directory")
    parser.add_argument("--outfile", type=str, required=True, help="path to the output plot")
    args = parser.parse_args()

    ## Read CAMI OPAL result file
    cami_opal_res_path = Path(args.cami_opal_res)
    if cami_opal_res_path.is_file():
        cami_opal_res = pd.read_csv(cami_opal_res_path, sep="\t", header=0)
    else:
        raise Exception("CAMI OPAL result file does not exist.")

    ## Write YACHT OPAL result
    yacht_opal_dir_path = Path(args.yacht_opal_dir)
    if yacht_opal_dir_path.is_dir():
        yacht_opal_res = pd.concat([pd.read_csv(file, sep='\t', header=0) for file in yacht_opal_dir_path.glob("*/results.tsv")])
    else:
        raise Exception("YACHT OPAL result directory does not exist.")
    

    ## Merge CAMI OPAL result and YAML OPAL result
    # remove some inappropriate tools
    if 'rhizosphere_data' in str(cami_opal_res_path):
        removed_index = cami_opal_res["tool"].str.contains('(nano)') | cami_opal_res["tool"].str.contains('(pb)') | cami_opal_res["tool"].str.contains('k21') | cami_opal_res["tool"].str.contains('k51') | cami_opal_res["tool"].str.contains('MetaPhlAn 2.9.21') | cami_opal_res["tool"].str.contains('mOTUs cami1')
        cami_opal_res = cami_opal_res.loc[~removed_index, :].reset_index(drop=True)
    elif 'marine_data' in str(cami_opal_res_path) or 'strain_madness_data' in str(cami_opal_res_path):
        cami_opal_res = cami_opal_res.loc[~cami_opal_res["tool"].isin(['mOTUv1', 'mOTUv2b', 'mOTUv2c', 'mOTUv2d', 'mOTUv2e', 'mOTUv2f', 'mOTUv2g', 'mOTUv2h', 'mOTUv2i']),:].reset_index(drop=True)
    else:
        raise Exception("The script now only supports the following datasets: rhizosphere_data, marine_data, strain_madness_data.")
    
    yacht_opal_res.loc[yacht_opal_res['tool'] != 'Gold standard', 'tool'] = 'YACHT k31'
    yacht_opal_res = yacht_opal_res.loc[yacht_opal_res['tool'] != 'Gold standard', :].reset_index(drop=True)
    combine_res = pd.concat([cami_opal_res, yacht_opal_res], ignore_index=True).reset_index(drop=True)
    
    # Specify the custom order
    temp_list = list(combine_res['tool'].unique())
    temp_list.remove('Gold standard')
    temp_list = sorted(temp_list, key=lambda x: x.upper()[0])
    order = ['Gold standard'] + temp_list

    # Convert the column to a categorical type with the specified order
    combine_res['tool'] = pd.Categorical(combine_res['tool'], categories=order, ordered=True)

    # Sort by the custom order
    combine_res = combine_res.sort_values('tool').reset_index(drop=True)
    
    # Get unique metrics
    metrics = combine_res['metric'].unique()
    # remove some metrics
    metrics = metrics[~pd.Series(metrics).isin(['Sum of abundances', 'Unweighted UniFrac error', 'Weighted UniFrac error', 'Unweighted UniFrac (CAMI)','Weighted UniFrac (CAMI)'])]

    # Create a combined plot
    fig, axes = plt.subplots(4, 3, figsize=(20, 15))
    fig.subplots_adjust(top=0.92)  # Adjust the top to make space for legend

    for idx, metric in enumerate(metrics):
        ax = axes[idx // 3, idx % 3]
        sns.boxplot(x="rank", y="value", hue="tool", data=combine_res[combine_res['metric'] == metric], order=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'], ax=ax)
        ax.set_title(metric)
        ax.set_xticklabels(['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'], rotation=20)
        ax.get_legend().remove()
        
    # Create a combined legend at the top
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=len(combine_res['tool'].unique()), bbox_to_anchor=(0.5, 1.02))

    plt.tight_layout()
    ## save figure to png
    plt.savefig(args.outfile, dpi=300, bbox_inches='tight')