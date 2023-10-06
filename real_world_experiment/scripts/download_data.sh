#! /bin/bash
# Download real world data from study (https://www.nature.com/articles/s41591-020-1105-z)

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: download_data.sh <yacht_reproducibles_dir> <cpu_num>"
    exit 1
fi
yacht_reproducibles_dir=$1
yacht_repo_loc=${yacht_reproducibles_dir}/YACHT
real_world_experiment_dir=${yacht_reproducibles_dir}/real_world_experiment
cpu_num=$2
data_dir=$real_world_experiment_dir/data

## define a bash function
function download_fastq() {
    dirname=`echo $2 | cut -d';' -f 1 | cut -d'/' -f 6`;
    # create the directory if not exist
    if [ ! -d $1/fastq/$dirname ]; then
        mkdir $1/fastq/$dirname;
    fi
    cd $1/fastq/$dirname;
    echo $2 | cut -d';' -f 1 | wget -i -;
    echo $2 | cut -d';' -f 2 | wget -i -;
}
export -f download_fastq

# download the data from https://www.ebi.ac.uk/ena/browser/view/PRJEB14847
echo "Downloading the sample data"
cd $data_dir
# ## extract the bioproject info
# if [ ! -f $data_dir/ncbi_bioproject_info.csv ]; then
#     esearch -db sra -query 'PRJNA558701' | efetch -format runinfo > ncbi_bioproject_info.csv
# fi
## download the fastq files
if [ ! -d $data_dir/fastq ]; then
    mkdir -p $data_dir/fastq;
    awk -F'\t' '$2=="Q" {print $1}' phaseIII_sample_list.tsv | awk -F'\t' 'NR==FNR {samples[$1]=1; next} {for (sample in samples) if ($7 ~ sample) print $8}' - filereport_read_run_PRJEB14847.tsv | parallel -j $cpu_num --link download_fastq $(pwd) {};
fi


