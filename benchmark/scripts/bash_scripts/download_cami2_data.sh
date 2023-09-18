#! /bin/bash
# Download CAMI2 data from https://data.cami-challenge.org/participate

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: download_cami2_data.sh <benchmark_dir> <cpu_num>"
    exit 1
fi
benchmark_dir=$1
cpu_num=$2

# create a CAMI_data directory for the output
if [ ! -d $benchmark_dir/CAMI_data ]; then
    mkdir $benchmark_dir/CAMI_data
fi

# download the data seperately for Rhizosphere challenge, Clinical pathogen detection challenge, Challenge Marine Dataset, Strain Madness Dataset
# Rhizosphere challenge
echo "Downloading Rhizosphere challenge data"
if [ ! -d $benchmark_dir/CAMI_data/rhizosphere_data ]; then
    mkdir $benchmark_dir/CAMI_data/rhizosphere_data;
    ## download the simulated data
    parallel -j $cpu_num wget -P $benchmark_dir/CAMI_data/rhizosphere_data https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/rhimgCAMI2_sample_{}.tar.gz ::: {0..20};
    ## Extract the fq.gz files
    cd $benchmark_dir/CAMI_data/rhizosphere_data;
    ## backup the raw files
    if [ ! -d ./backup ]; then
        mkdir ./backup;
    fi
    mv ./*.tar.gz ./backup;
    ## extract the fq.gz files
    parallel -j $cpu_num tar zxvf ./backup/rhimgCAMI2_sample_{}_reads.tar.gz ::: {0..20};
    parallel -j $cpu_num mv ./simulation_short_read/2019.09.27_13.59.10_sample_{}/reads/anonymous_reads.fq.gz ./rhimgCAMI2_sample_{}.fq.gz ::: {0..20};
    rm -r ./simulation_short_read;
    ## download the OPAL result
    wget -P $benchmark_dir/CAMI_data/rhizosphere_data/opal_results https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/rhizosphere_dataset/results/OPAL_short_long_noplasmids_normalized_filtered/results.tsv;
    ## download the ground truth profile
    wget -P $benchmark_dir/CAMI_data/rhizosphere_data/ground_truth https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/rhizosphere_dataset/data/ground_truth/gs_rhizosphere.filtered.profile;
    # Split the ground truth profile into individual files
    cd $benchmark_dir/CAMI_data/rhizosphere_data/ground_truth
    less gs_rhizosphere.filtered.profile | sed 's/short_read_//' | awk '/@SampleID:/{filename=substr($0, length("@SampleID:") + 1) ".profile"} {print >filename}'
        
    ## download raw genome files and their metadata
    wget -P $benchmark_dir/CAMI_data/rhizosphere_data/genomes https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/rhimgCAMI2_setup.tar.gz;
    tar zxvf $benchmark_dir/CAMI_data/rhizosphere_data/genomes/rhimgCAMI2_setup.tar.gz -C $benchmark_dir/CAMI_data/rhizosphere_data/genomes;
    mv $benchmark_dir/CAMI_data/rhizosphere_data/genomes/simulation_short_read/metadata.tsv $benchmark_dir/CAMI_data/rhizosphere_data/genomes;
    mv $benchmark_dir/CAMI_data/rhizosphere_data/genomes/simulation_short_read/genome_to_id.tsv $benchmark_dir/CAMI_data/rhizosphere_data/genomes;
    rm -r $benchmark_dir/CAMI_data/rhizosphere_data/genomes/simulation_short_read  $benchmark_dir/CAMI_data/rhizosphere_data/genomes/rhimgCAMI2_setup.tar.gz;
    wget -P $benchmark_dir/CAMI_data/rhizosphere_data/genomes https://frl.publisso.de/data/frl:6425521/plant_associated/rhimgCAMI2_genomes.tar.gz;
    tar zxvf $benchmark_dir/CAMI_data/rhizosphere_data/genomes/rhimgCAMI2_genomes.tar.gz -C $benchmark_dir/CAMI_data/rhizosphere_data/genomes;
    mv $benchmark_dir/CAMI_data/rhizosphere_data/genomes/source_genomes/*.fasta $benchmark_dir/CAMI_data/rhizosphere_data/genomes;
    rm -r $benchmark_dir/CAMI_data/rhizosphere_data/genomes/source_genomes $benchmark_dir/CAMI_data/rhizosphere_data/genomes/rhimgCAMI2_genomes.tar.gz;
    ## combine metadata
    python $benchmark_dir/scripts/python_scripts/combine_metadata_files.py --metadata $benchmark_dir/CAMI_data/rhizosphere_data/genomes/metadata.tsv --genome_to_id $benchmark_dir/CAMI_data/rhizosphere_data/genomes/genome_to_id.tsv --ground_truth_dir $benchmark_dir/CAMI_data/rhizosphere_data/ground_truth;
fi


# Clinical pathogen detection challenge
echo "Downloading Clinical pathogen detection challenge data"
if [ ! -d $benchmark_dir/CAMI_data/pathogen_detection_data ]; then
    mkdir $benchmark_dir/CAMI_data/pathogen_detection_data;
    wget -P $benchmark_dir/CAMI_data/pathogen_detection_data https://frl.publisso.de/data/frl:6425521/patmgCAMI2.tar.gz;
    tar zxvf $benchmark_dir/CAMI_data/pathogen_detection_data/patmgCAMI2.tar.gz -C $benchmark_dir/CAMI_data/pathogen_detection_data;
    mv $benchmark_dir/CAMI_data/pathogen_detection_data/patmg_CAMI2/* $benchmark_dir/CAMI_data/pathogen_detection_data;
    rm -r $benchmark_dir/CAMI_data/pathogen_detection_data/patmg_CAMI2;
    rm $benchmark_dir/CAMI_data/pathogen_detection_data/patmgCAMI2.tar.gz;
fi


# Challenge Marine Dataset
echo "Downloading Challenge Marine Dataset data"
if [ ! -d $benchmark_dir/CAMI_data/marine_data ]; then
    mkdir $benchmark_dir/CAMI_data/marine_data
    parallel -j $cpu_num wget -P $benchmark_dir/CAMI_data/marine_data https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_sample_{}_reads.tar.gz ::: {0..9};
    ## Extract the fq.gz files
    cd $benchmark_dir/CAMI_data/marine_data;
    ## backup the raw files
    if [ ! -d ./backup ]; then
        mkdir ./backup;
    fi
    mv ./*.tar.gz ./backup;
    ## extract the fq.gz files
    parallel -j $cpu_num tar zxvf ./backup/marmgCAMI2_sample_{}_reads.tar.gz ::: {0..9};
    parallel -j $cpu_num mv ./simulation_short_read/2018.08.15_09.49.32_sample_{}/reads/anonymous_reads.fq.gz ./marmgCAMI2_sample_{}.fq.gz ::: {0..9};
    rm -r ./simulation_short_read;
    # download the OPAL result
    wget -P $benchmark_dir/CAMI_data/marine_data/opal_results https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/marine_dataset/results/OPAL_short_long_noplasmids_normalized_filtered/results.tsv;
    # download the ground truth profile
    wget -P $benchmark_dir/CAMI_data/marine_data/ground_truth https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/marine_dataset/data/ground_truth/gs_marine_short.filtered.profile;
    # Split the ground truth profile into individual files
    cd $benchmark_dir/CAMI_data/marine_data/ground_truth
    less gs_marine_short.profile | sed 's/short_read_//' | awk '/@SampleID:/{filename=substr($0, length("@SampleID:") + 1) ".profile"} {print >filename}'
        
    ## download raw genome files and their metadata
    wget -P $benchmark_dir/CAMI_data/marine_data/genomes https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_setup.tar.gz;
    tar zxvf $benchmark_dir/CAMI_data/marine_data/genomes/marmgCAMI2_setup.tar.gz -C $benchmark_dir/CAMI_data/marine_data/genomes;
    mv $benchmark_dir/CAMI_data/marine_data/genomes/simulation_short_read/metadata.tsv $benchmark_dir/CAMI_data/marine_data/genomes;
    mv $benchmark_dir/CAMI_data/marine_data/genomes/simulation_short_read/genome_to_id.tsv $benchmark_dir/CAMI_data/marine_data/genomes;
    rm -r $benchmark_dir/CAMI_data/marine_data/genomes/simulation_short_read  $benchmark_dir/CAMI_data/marine_data/genomes/marmgCAMI2_setup.tar.gz;
    wget -P $benchmark_dir/CAMI_data/marine_data/genomes https://frl.publisso.de/data/frl:6425521/marine/marmgCAMI2_genomes.tar.gz;
    tar zxvf $benchmark_dir/CAMI_data/marine_data/genomes/rhimgCAMI2_genomes.tar.gz -C $benchmark_dir/CAMI_data/marine_data/genomes;
    mv $benchmark_dir/CAMI_data/marine_data/genomes/simulation_short_read/genomes/*.fasta $benchmark_dir/CAMI_data/marine_data/genomes;
    rm -r $benchmark_dir/CAMI_data/marine_data/genomes/simulation_short_read $benchmark_dir/CAMI_data/marine_data/genomes/marmgCAMI2_genomes.tar.gz;
    ## combine metadata
    python $benchmark_dir/scripts/python_scripts/combine_metadata_files.py --metadata $benchmark_dir/CAMI_data/marine_data/genomes/metadata.tsv --genome_to_id $benchmark_dir/CAMI_data/marine_data/genomes/genome_to_id.tsv --ground_truth_dir $benchmark_dir/CAMI_data/marine_data/ground_truth;
fi

# Strain Madness Dataset
echo "Downloading Strain Madness Dataset data"
if [ ! -d $benchmark_dir/CAMI_data/strain_madness_data ]; then
    mkdir $benchmark_dir/CAMI_data/strain_madness_data;
    parallel -j $cpu_num wget -P $benchmark_dir/CAMI_data/strain_madness_data https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_sample_{}_reads.tar.gz ::: {0..99};
    ## Extract the fq.gz files
    cd $benchmark_dir/CAMI_data/strain_madness_data;
    ## backup the raw files
    if [ ! -d ./backup ]; then
        mkdir ./backup;
    fi
    mv ./*.tar.gz ./backup;
    ## extract the fq.gz files
    parallel -j $cpu_num tar zxvf ./backup/strmgCAMI2_sample_{}_reads.tar.gz ::: {0..99};
    parallel -j $cpu_num mv ./short_read/2018.09.07_11.43.52_sample_{}/reads/anonymous_reads.fq.gz ./strmgCAMI2_sample_{}.fq.gz ::: {0..99};
    rm -r ./short_read;
    # download the OPAL result
    wget -P $benchmark_dir/CAMI_data/strain_madness_data/opal_results https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/strain_madness_dataset/results/OPAL_default_short_read_samples/results.tsv;
    # download the ground truth profile
    wget -P $benchmark_dir/CAMI_data/strain_madness_data/ground_truth https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/strain_madness_dataset/data/ground_truth/gs_strain_madness_short_long.profile;
    # Split the ground truth profile into individual files
    cd $benchmark_dir/CAMI_data/strain_madness_data/ground_truth
    less gs_strain_madness_short_long.profile | sed 's/short_read_//' | awk '/@SampleID:/{filename=substr($0, length("@SampleID:") + 1) ".profile"} {print >filename}'


    ## download raw genome files and their metadata
    wget -P $benchmark_dir/CAMI_data/strain_madness_data/genomes https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_setup.tar.gz;
    tar zxvf $benchmark_dir/CAMI_data/strain_madness_data/genomes/strmgCAMI2_setup.tar.gz -C $benchmark_dir/CAMI_data/strain_madness_data/genomes;
    mv $benchmark_dir/CAMI_data/strain_madness_data/genomes/short_read/metadata.tsv $benchmark_dir/CAMI_data/strain_madness_data/genomes;
    mv $benchmark_dir/CAMI_data/strain_madness_data/genomes/short_read/genome_to_id.tsv $benchmark_dir/CAMI_data/strain_madness_data/genomes;
    rm -r $benchmark_dir/CAMI_data/strain_madness_data/genomes/short_read  $benchmark_dir/CAMI_data/strain_madness_data/genomes/strmgCAMI2_setup.tar.gz;
    wget -P $benchmark_dir/CAMI_data/strain_madness_data/genomes https://frl.publisso.de/data/frl:6425521/strain/strmgCAMI2_genomes.tar.gz;
    tar zxvf $benchmark_dir/CAMI_data/strain_madness_data/genomes/strmgCAMI2_genomes.tar.gz -C $benchmark_dir/CAMI_data/strain_madness_data/genomes;
    mv $benchmark_dir/CAMI_data/strain_madness_data/genomes/short_read/source_genomes/*.fasta $benchmark_dir/CAMI_data/strain_madness_data/genomes;
    rm -r $benchmark_dir/CAMI_data/strain_madness_data/genomes/short_read $benchmark_dir/CAMI_data/strain_madness_data/genomes/strmgCAMI2_genomes.tar.gz;
    ## combine metadata
    python $benchmark_dir/scripts/python_scripts/combine_metadata_files.py --metadata $benchmark_dir/CAMI_data/strain_madness_data/genomes/metadata.tsv --genome_to_id $benchmark_dir/CAMI_data/strain_madness_data/genomes/genome_to_id.tsv --ground_truth_dir $benchmark_dir/CAMI_data/strain_madness_data/ground_truth;
fi