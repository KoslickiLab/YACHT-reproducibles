"""
This script is used to generate the plots for comparing YACHT with other SOTA tools.
"""


import os
from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.markers as mmarkers
import matplotlib.patches as patches
import seaborn as sns
import math

def get_all_markers():
    return list(mmarkers.MarkerStyle.markers.keys())

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
        yacht_opal_res_dict = {}
        for file in yacht_opal_dir_path.glob("*/results.tsv"):
            key = str(file).split('/')[-2].split('_')[0]
            if key not in yacht_opal_res_dict:
                yacht_opal_res_dict[key] = []
            yacht_opal_res_dict[key] += [pd.read_csv(file, sep='\t', header=0)]
        for key in yacht_opal_res_dict:
            yacht_opal_res_dict[key] = pd.concat(yacht_opal_res_dict[key], ignore_index=True)
    else:
        raise Exception("YACHT OPAL result directory does not exist.")
    

    ## Merge CAMI OPAL result and YAML OPAL result
    # remove some inappropriate tools
    if 'rhizosphere_data' in str(cami_opal_res_path):
        pattern = '|'.join(['(nano)', '(pb)', 'k21', 'k51', 'MetaPhlAn 2.9.21', 'mOTUs cami1'])
        removed_index = cami_opal_res["tool"].str.contains(pattern)
        cami_opal_res = cami_opal_res.loc[~removed_index, :].reset_index(drop=True)
        cami_opal_res['tool'] = cami_opal_res['tool'].str.replace(' (sr)', '')
    elif 'marine_data' in str(cami_opal_res_path):
        cami_opal_res = cami_opal_res.loc[~cami_opal_res["tool"].isin(['Metalign 0.6.2 avg', 'LSHVec illumina', 'LSHVec pacbio', 'mOTUs 2.0.1_1', 'mOTUs 2.5.1_11', 
                                                                       'mOTUs 2.5.1_3','mOTUs 2.5.1_4', 'mOTUs 2.5.1_5', 'mOTUs 2.5.1_6', 'mOTUs 2.5.1_7', 
                                                                       'mOTUs 2.5.1_8', 'mOTUs 2.5.1_9', 'mOTUs 2.5.1_10',
                                                                       'DUDes cami1', 'FOCUS cami1', 'MetaPhlAn cami1', 'mOTUs cami1', 'TIPP cami1']),:].reset_index(drop=True)
        cami_opal_res['tool'] = cami_opal_res['tool'].str.replace('LSHVec gsa', 'LSHVec')
        cami_opal_res['tool'] = cami_opal_res['tool'].str.replace('mOTUs 2.5.1_2', 'mOTUs 2.5.1')
    elif 'strain_madness_data' in str(cami_opal_res_path):
        cami_opal_res = cami_opal_res.loc[~cami_opal_res["tool"].isin(['mOTUv1', 'mOTUv2b', 'mOTUv2c', 'mOTUv2d', 
                                                                       'mOTUv2e', 'mOTUv2f', 'mOTUv2g', 'mOTUv2h', 'mOTUv2i']),:].reset_index(drop=True)
        cami_opal_res['tool'] = cami_opal_res['tool'].str.replace('mOTUv2a', 'mOTUv2')
        cami_opal_res['tool'] = cami_opal_res['tool'].str.replace('Metln', 'Metalign')
        cami_opal_res['tool'] = cami_opal_res['tool'].str.replace('MP', 'MetaPhlAn')
    else:
        raise Exception("The script now only supports the following datasets: rhizosphere_data, marine_data, strain_madness_data.")
    
    for key in yacht_opal_res_dict:
        temp = yacht_opal_res_dict[key]
        temp.loc[temp['tool'] != 'Gold standard', 'tool'] = f'YACHT k31 (min_{key})'
        temp = temp.loc[temp['tool'] != 'Gold standard', :].reset_index(drop=True)
        yacht_opal_res_dict[key] = temp
    yacht_opal_res = pd.concat(yacht_opal_res_dict.values(), ignore_index=True).reset_index(drop=True)
    combine_res = pd.concat([cami_opal_res, yacht_opal_res], ignore_index=True).reset_index(drop=True)
    sebuset_combine_res = combine_res.loc[(combine_res['tool'] != 'Gold standard') & (combine_res['metric'].isin(['Completeness', 'Purity', 'F1 score'])) & (combine_res['rank'] != 'strain'),:].reset_index(drop=True)
    
    # Specify the custom order
    temp_list = list(combine_res['tool'].unique())
    temp_list.remove('Gold standard')
    order1 = sorted([x for x in temp_list if 'YACHT' not in x], key=lambda x: x.upper()[0]) + \
        sorted([x for x in temp_list if 'YACHT' in x], key=lambda x: float(x.replace('YACHT k31 (min_coverage', '').replace(')', '')), reverse=True)
    order2 = ['Completeness', 'Purity', 'F1 score']

    # Convert the column to a categorical type with the specified order
    sebuset_combine_res['tool'] = pd.Categorical(sebuset_combine_res['tool'], categories=order1, ordered=True)
    sebuset_combine_res['metric'] = pd.Categorical(sebuset_combine_res['metric'], categories=order2, ordered=True)

    # Sort by the custom order
    sebuset_combine_res = sebuset_combine_res.sort_values(['tool','metric']).reset_index(drop=True)
    
    # Get unique tools
    tools = sebuset_combine_res['tool'].unique()

    # # Create a combined plot
    # fig, axes = plt.subplots(math.ceil(len(sebuset_combine_res['tool'].unique()) / 3), 3, figsize=(15, math.ceil(len(tools)/3)*3))
    # fig.subplots_adjust(top=0.98)  # Adjust the top to make space for legend

    # for idx, tool in enumerate(tools):
    #     ax = axes[idx // 3, idx % 3]
    #     sns.lineplot(x="rank", y="value", hue="metric", errorbar='sd', data=sebuset_combine_res[sebuset_combine_res['tool'] == tool], ax=ax)
    #     ax.set_title(tool)
    #     ax.set_xlabel("")
    #     ax.set_xticklabels(['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'], rotation=30)
    #     ax.set_ylim(0, 1.1)
    #     ax.get_legend().remove()

    # ## remove the empty subfigures
    # for i in range(len(tools), math.ceil(len(sebuset_combine_res['tool'].unique()) / 3) * 3):
    #     fig.delaxes(axes[i // 3, i % 3])

    # # Create a combined legend at the top
    # handles, labels = ax.get_legend_handles_labels()
    # fig.legend(handles, labels, loc='upper center', ncol=len(sebuset_combine_res['metric'].unique()), bbox_to_anchor=(0.5, 1.01))

    # # Add the title above the legend
    # if args.fig_title is not None:
    #     fig.suptitle(args.fig_title, fontsize=12, y=1.04)

    # plt.tight_layout()
    # ## save figure
    # plt.savefig(os.path.join(args.outdir, args.filename+'_lineplot.png'), dpi=300, bbox_inches='tight')

    # draw a box plot with completeness, purity and f1 score
    for metric in ['Completeness', 'Purity', 'F1 score']:
        fig, ax = plt.subplots(figsize=(10, 6))

        # Define box properties to have no facecolor
        boxprops = {'facecolor': 'none'}

        sns.boxplot(x='tool', y='value', data=sebuset_combine_res.query(f'rank == "species" and metric == "{metric}"').reset_index(drop=True), boxprops=boxprops, ax=ax)

        # Hight the yacht with grey
        boxes = [patch for patch in ax.get_children() if isinstance(patch, patches.PathPatch)]
        for i in range(1, 6):
            boxes[-i].set_facecolor('grey')

        ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
        ax.set_xlabel('')
        ax.set_ylabel(metric)
        ax.set_ylim(0, 1.1)
        if args.fig_title is not None:
            ax.set_title(args.fig_title, fontsize=12)
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        ## save figure
        plt.savefig(os.path.join(args.outdir, args.filename+f"_{metric.replace(' ','').lower()}_boxplot.png"), dpi=300, bbox_inches='tight')

    ## draw a scatter plot with completeness vs purity and error bar
    sebuset_combine_res = sebuset_combine_res.loc[sebuset_combine_res['metric'] != 'F1 score', :].reset_index(drop=True)
    # Group by tool, rank and metric
    grouped = sebuset_combine_res.groupby(['tool', 'rank', 'metric'])['value'].agg(['mean', 'std'])

    # Define different markers for tools
    markers = ['o', 'v', '^', '<', '>', 's', 'D', '*', 'p', 'H', 'X'] # get_all_markers()
    markers += markers[::-1]

    # Define different colors
    nonyacht_colors = [x for index, x in enumerate(sns.color_palette("Set3", 11)) if index not in [1, 8]]
    nonyacht_colors += nonyacht_colors[1:] + [nonyacht_colors[0]]
    yacht_color = 'grey'

    ranks = sebuset_combine_res['rank'].unique()
    num_ranks = len(ranks)

    # Create a single scatter plot
    fig, ax = plt.subplots(figsize=(8, 5))

    rank = 'species'
    temp_tool_list = sebuset_combine_res['tool'].unique()
    non_yacht_tools = [x for x in temp_tool_list if 'YACHT' not in x]
    yacht_tools = [x for x in temp_tool_list if 'YACHT' in x]
    
    for j, tool in enumerate(non_yacht_tools):
        x = grouped.loc[tool, rank].xs('Purity')['mean']
        xerr = grouped.loc[tool, rank].xs('Purity')['std']
        y = grouped.loc[tool, rank].xs('Completeness')['mean']
        yerr = grouped.loc[tool, rank].xs('Completeness')['std']
        
        line = ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=markers[j], label=tool, ecolor=nonyacht_colors[j], color=nonyacht_colors[j], capsize=5)
    
    for i, tool in enumerate(yacht_tools):
        i = i + len(non_yacht_tools)
        x = grouped.loc[tool, rank].xs('Purity')['mean']
        xerr = grouped.loc[tool, rank].xs('Purity')['std']
        y = grouped.loc[tool, rank].xs('Completeness')['mean']
        yerr = grouped.loc[tool, rank].xs('Completeness')['std']
        
        line = ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=markers[i], label=tool, ecolor=yacht_color, color=yacht_color, capsize=5)
    
    ax.set_xlabel('Purity')
    ax.set_ylabel('Completeness')
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 1.1)
    if args.fig_title is not None:
        ax.set_title(args.fig_title, fontsize=12)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    # # Create a combined plot
    # fig, axes = plt.subplots(math.ceil(num_ranks / 4), 4, figsize=(15, 8))
    # fig.subplots_adjust(top=0.98)  # Adjust the top to make space for legend

    # for idx, rank in enumerate(ranks):
    #     ax = axes[idx // 4, idx % 4]
    #     for j, tool in enumerate(sebuset_combine_res['tool'].unique()):
    #         x = grouped.loc[tool, rank].xs('Purity')['mean']
    #         xerr = grouped.loc[tool, rank].xs('Purity')['std']
    #         y = grouped.loc[tool, rank].xs('Completeness')['mean']
    #         yerr = grouped.loc[tool, rank].xs('Completeness')['std']
            
    #         line = ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=markers[j], label=tool, capsize=5)
        
    #     ax.set_xlabel('Purity')
    #     ax.set_ylabel('Completeness')
    #     ax.set_xlim(0, 1.1)
    #     ax.set_ylim(0, 1.1)
    #     ax.set_title(f'{rank}')

    # ## remove the empty subfigures
    # for i in range(len(ranks), math.ceil(num_ranks / 4) * 4):
    #     fig.delaxes(axes[i // 4, i % 4])

    # # Create a combined legend at the top
    # handles, labels = ax.get_legend_handles_labels()
    # fig.legend(handles, labels, loc='upper center', ncol=len(sebuset_combine_res['tool'].unique()), bbox_to_anchor=(0.5, 1.01))

    # # Add the title above the legend
    # if args.fig_title is not None:
    #     fig.suptitle(args.fig_title, fontsize=12, y=1.04)
        
    plt.tight_layout()
    ## save figure
    plt.savefig(os.path.join(args.outdir, args.filename+'_scatterplot.png'), dpi=300, bbox_inches='tight')
