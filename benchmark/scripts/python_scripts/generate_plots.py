"""
This script is used to generate the plots for comparing YACHT with other SOTA tools.
"""


import os
from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.markers as mmarkers
import seaborn as sns
import math
import pickle

def get_all_markers():
    return list(mmarkers.MarkerStyle.markers.keys())

# sourcery skip: raise-specific-error
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is used to generate the plots for comparing YACHT with other SOTA tools.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cami_opal_res", type=str, required=True, help="path to the CAMI reported OPAL result file")
    parser.add_argument("--yacht_opal_dir", type=str, required=True, help="path to the YACHT OPAL result directory")
    parser.add_argument("--outdir", type=str, required=True, help="path to the output directory")
    parser.add_argument("--filename", type=str, required=False, help="output filename without suffix", default="result")
    parser.add_argument("--fig_title", type=str, required=False, help="figure title", default=None)
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
    sebuset_combine_res = combine_res.loc[(combine_res['tool'] != 'Gold standard') & (combine_res['metric'].isin(['Completeness', 'Purity', 'F1 score'])) & (combine_res['rank'] != 'strain'),:].reset_index(drop=True)
    
    # Specify the custom order
    temp_list = list(combine_res['tool'].unique())
    temp_list.remove('Gold standard')
    order1 = sorted(temp_list, key=lambda x: x.upper()[0])
    order2 = ['Completeness', 'Purity', 'F1 score']

    # Convert the column to a categorical type with the specified order
    sebuset_combine_res['tool'] = pd.Categorical(sebuset_combine_res['tool'], categories=order1, ordered=True)
    sebuset_combine_res['metric'] = pd.Categorical(sebuset_combine_res['metric'], categories=order2, ordered=True)

    # Sort by the custom order
    sebuset_combine_res = sebuset_combine_res.sort_values(['tool','metric']).reset_index(drop=True)
    
    # Get unique tools
    tools = sebuset_combine_res['tool'].unique()

    # Create a combined plot
    fig, axes = plt.subplots(math.ceil(len(sebuset_combine_res['tool'].unique()) / 3), 3, figsize=(15, 8))
    fig.subplots_adjust(top=0.98)  # Adjust the top to make space for legend

    for idx, tool in enumerate(tools):
        ax = axes[idx // 3, idx % 3]
        sns.lineplot(x="rank", y="value", hue="metric", data=sebuset_combine_res[sebuset_combine_res['tool'] == tool], ax=ax)
        ax.set_title(tool)
        ax.set_xlabel("")
        ax.set_xticklabels(['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'], rotation=30)
        ax.get_legend().remove()

    ## remove the empty subfigures
    for i in range(len(tools), math.ceil(len(sebuset_combine_res['tool'].unique()) / 3) * 3):
        fig.delaxes(axes[i // 3, i % 3])

    # Create a combined legend at the top
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=len(sebuset_combine_res['metric'].unique()), bbox_to_anchor=(0.5, 1.01))

    # Add the title above the legend
    if args.fig_title is not None:
        fig.suptitle(args.fig_title, fontsize=12, y=1.04)

    plt.tight_layout()
    ## save figure
    plt.savefig(os.path.join(args.outdir, args.filename+'_lineplot.png'), dpi=300, bbox_inches='tight')

    ## draw a scatter plot with completeness vs purity and error bar
    sebuset_combine_res = sebuset_combine_res.loc[sebuset_combine_res['metric'] != 'F1 score', :].reset_index(drop=True)
    # Group by tool, rank and metric
    grouped = sebuset_combine_res.groupby(['tool', 'rank', 'metric'])['value'].agg(['mean', 'std'])

    # Define different markers for tools
    markers = ['o', 's', 'D', '^', 'v', '*', 'p', 'H'] # get_all_markers()

    ranks = sebuset_combine_res['rank'].unique()
    num_ranks = len(ranks)

    # Create a combined plot
    fig, axes = plt.subplots(math.ceil(num_ranks / 4), 4, figsize=(15, 8))
    fig.subplots_adjust(top=0.98)  # Adjust the top to make space for legend

    for idx, rank in enumerate(ranks):
        ax = axes[idx // 4, idx % 4]
        for j, tool in enumerate(sebuset_combine_res['tool'].unique()):
            x = grouped.loc[tool, rank].xs('Purity')['mean']
            xerr = grouped.loc[tool, rank].xs('Purity')['std']
            y = grouped.loc[tool, rank].xs('Completeness')['mean']
            yerr = grouped.loc[tool, rank].xs('Completeness')['std']
            
            line = ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=markers[j], label=tool, capsize=5)
        
        ax.set_xlabel('Purity')
        ax.set_ylabel('Completeness')
        ax.set_title(f'{rank}')
        

    ## remove the empty subfigures
    for i in range(len(ranks), math.ceil(num_ranks / 4) * 4):
        fig.delaxes(axes[i // 4, i % 4])

    # Create a combined legend at the top
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=len(sebuset_combine_res['tool'].unique()), bbox_to_anchor=(0.5, 1.01))

    # Add the title above the legend
    if args.fig_title is not None:
        fig.suptitle(args.fig_title, fontsize=12, y=1.04)
        
    plt.tight_layout()
    ## save figure
    plt.savefig(os.path.join(args.outdir, args.filename+'_scatterplot.png'), dpi=300, bbox_inches='tight')
