#! /bin/bash
# run YACHT algorithm on a real-world data (https://www.ebi.ac.uk/ena/browser/view/PRJEB39681)
## Please be sure you have downloaded the data (e.g., falstq files and metadata) and put them in ~/YACHT-reproducibles/real_world_experiment/data

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run_YACHT.sh <yacht_reproducibles_dir> <cpu_num>"
    exit 1
fi
yacht_reproducibles_dir=$1
yacht_repo_loc=${yacht_reproducibles_dir}/YACHT
real_world_experiment_dir=${yacht_reproducibles_dir}/real_world_experiment
cpu_num=$2

## define a bash function
function download_microbe_assembly() {
    cd $1/pathogen_detection_reference/microbe_genomes;
    if [ ! -e $1/pathogen_detection_reference/microbe_genomes/${2}.fna ]; then
        mkdir $1/pathogen_detection_reference/microbe_genomes/$2;
        cd $1/pathogen_detection_reference/microbe_genomes/$2;
        datasets download genome taxon $2 --assembly-level complete --assembly-source all --no-progressbar --tax-exact-match --filename ${2}.zip;
        if [ ! -e $1/pathogen_detection_reference/microbe_genomes/$2/${2}.zip ]; then
            datasets download genome taxon $2 --assembly-source all --no-progressbar --tax-exact-match --filename ${2}.zip;
            if [ -e $1/pathogen_detection_reference/microbe_genomes/$2/${2}.zip ]; then
                unzip ${2}.zip;
                filename=`find ./ncbi_dataset -type f -name "*.fna" | sort -t'/' -k4 -r | head -1`
            else
                filename="";
            fi
        else
            unzip ${2}.zip;
            filename=`find ./ncbi_dataset -type f -name "*.fna" | sort -t'/' -k4 -r | head -1`
        fi
        if [ -n "$filename" ]; then
            echo "$2" >> $1/pathogen_detection_reference/microbe_check_list.txt;
            mv $filename $1/pathogen_detection_reference/microbe_genomes/${2}.fna;
        fi
        rm -rf $1/pathogen_detection_reference/microbe_genomes/$2;
    fi
}
export -f download_microbe_assembly

function download_fungi_assembly() {
    cd $1/pathogen_detection_reference/fungi_genomes;
    if [ ! -e $1/pathogen_detection_reference/fungi_genomes/${2}.fna ]; then
        mkdir $1/pathogen_detection_reference/fungi_genomes/$2;
        cd $1/pathogen_detection_reference/fungi_genomes/$2;
        datasets download genome taxon $2 --assembly-level complete --assembly-source all --no-progressbar --tax-exact-match --filename ${2}.zip;
        if [ ! -e $1/pathogen_detection_reference/fungi_genomes/$2/${2}.zip ]; then
            datasets download genome taxon $2 --assembly-source all --no-progressbar --tax-exact-match --filename ${2}.zip;
            if [ -e $1/pathogen_detection_reference/fungi_genomes/$2/${2}.zip ]; then
                unzip ${2}.zip;
                filename=`find ./ncbi_dataset -type f -name "*.fna" | sort -t'/' -k4 -r | head -1`
            else
                filename="";
            fi
        else
            unzip ${2}.zip;
            filename=`find ./ncbi_dataset -type f -name "*.fna" | sort -t'/' -k4 -r | head -1`
        fi
        if [ -n "$filename" ]; then
            echo "$2" >> $1/pathogen_detection_reference/fungi_check_list.txt;
            mv $filename $1/pathogen_detection_reference/fungi_genomes/${2}.fna;
        fi
        rm -rf $1/pathogen_detection_reference/fungi_genomes/$2;
    fi
}
export -f download_fungi_assembly

## define a bash function
function download_virus_assembly() {
    cd $1/pathogen_detection_reference/virus_genomes;
    if [ ! -d $1/pathogen_detection_reference/virus_genomes/$2 ]; then
        mkdir $1/pathogen_detection_reference/virus_genomes/$2;
        cd $1/pathogen_detection_reference/virus_genomes/$2;
        ncbi-genome-download -s refseq -F "fasta" -l "complete" -t $2 viral;
        if [ ! -d $1/pathogen_detection_reference/virus_genomes/$2/refseq ]; then
            ncbi-genome-download -s genbank -F "fasta" -l "complete" -t $2 viral;
            if [ -d $1/pathogen_detection_reference/virus_genomes/$2/genbank ]; then
                filename=`find ./genbank -type f -name "*.fna.gz" | head -1`
            else
                filename=""
            fi
        else
            filename=`find ./refseq -type f -name "*.fna.gz" | head -1`
        fi
        if [ -n "$filename" ]; then
            echo "$2" >> $1/pathogen_detection_reference/virus_check_list.txt;
            mv $filename $1/pathogen_detection_reference/virus_genomes/${2}.fna.gz;
        fi
        rm -rf $1/pathogen_detection_reference/virus_genomes/$2;
    else
        cd $1/pathogen_detection_reference/virus_genomes/$2;
        ncbi-genome-download -s refseq -F "fasta" -l "complete" -t $2 viral;
        if [ ! -d $1/pathogen_detection_reference/virus_genomes/$2/refseq ]; then
            ncbi-genome-download -s genbank -F "fasta" -l "complete" -t $2 viral;
            if [ -d $1/pathogen_detection_reference/virus_genomes/$2/genbank ]; then
                filename=`find ./genbank -type f -name "*.fna.gz" | head -1`
            else
                filename=""
            fi
        else
            filename=`find ./refseq -type f -name "*.fna.gz" | head -1`
        fi
        if [ -n "$filename" ]; then
            echo "$2" >> $1/pathogen_detection_reference/virus_check_list.txt;
            mv $filename $1/pathogen_detection_reference/virus_genomes/${2}.fna.gz;
            gunzip $1/pathogen_detection_reference/virus_genomes/${2}.fna.gz;
        fi
        rm -rf $1/pathogen_detection_reference/virus_genomes/$2;
    fi
}
export -f download_virus_assembly


# Build a pathogen database (data collected from BVBRC) used for CAMI2 Challenge
echo "Building a pathogen database (data collected from BVBRC) used for CAMI2 Challenge"
cd $yacht_repo_loc
if [ ! -d $yacht_repo_loc/pathogen_detection_reference ]; then
    mkdir $yacht_repo_loc/pathogen_detection_reference;
    cd $yacht_repo_loc/pathogen_detection_reference;
    # download non-viral genomes from BVBRC
    wget -O $yacht_repo_loc/pathogen_detection_reference/microbe_genome_metadata.tsv ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata;
    awk 'BEGIN{IGNORECASE=1;FS="\t";OFS="\t"} {print $1,$2,$4,$5,$18,$19,$20,$46,$64}' microbe_genome_metadata.tsv | cut -f 3 | sort -u > microbe_taxids.txt
    taxonkit lineage -c -r -L microbe_taxids.txt | awk -F'\t' '$3 ~ /species|strain/ {print $0}' > microbe_taxids_filtered.txt;
    if [ ! -d $yacht_repo_loc/pathogen_detection_reference/microbe_genomes ]; then
        mkdir $yacht_repo_loc/pathogen_detection_reference/microbe_genomes;
        cd $yacht_repo_loc/pathogen_detection_reference/microbe_genomes;
        less $yacht_repo_loc/pathogen_detection_reference/microbe_taxids_filtered.txt | cut -f 2 | parallel -j $cpu_num "datasets download genome taxon {} --assembly-level complete --assembly-source all --no-progressbar --tax-exact-match --filename '{}.zip'";
        for file in *.zip; do
            unzip $file;
            filename=`find ./ncbi_dataset -type f -name "*.fna" | sort -t'/' -k4 -r | head -1`
            taxid=`echo $file | sed 's/\.zip//'`
            if [ -n "$filename" ]; then
                echo "$taxid" >> $yacht_repo_loc/pathogen_detection_reference/microbe_check_list.txt;
                mv $filename ${taxid}.fna;    
            fi
            rm -rf $file ncbi_dataset README.md;
        done
        echo 1260 36470 120793 1344959 77643 480036 671232 502800 | tr ' ' '\n' > $yacht_repo_loc/pathogen_detection_reference/microbe_taxid_add_list.txt
        less $yacht_repo_loc/pathogen_detection_reference/microbe_taxid_add_list.txt | parallel -j $cpu_num --link download_microbe_assembly $yacht_repo_loc {};
    fi

    cd $yacht_repo_loc/pathogen_detection_reference;
    # download viral genomes from BVBRC
    wget -O $yacht_repo_loc/pathogen_detection_reference/virus_genome_metadata.tsv ftp://ftp.bvbrc.org/viruses/genome_metadata;
    awk 'BEGIN{IGNORECASE=1;FS="\t";OFS="\t"}{print $1,$2,$4,$5,$18,$19,$20,$46,$64}' virus_genome_metadata.tsv | cut -f 3 | sort -u > virus_taxids.txt
    taxonkit lineage -c -r -L virus_taxids.txt | awk -F'\t' '$3 ~ /rank|species|strain/ {print $0}' > virus_taxids_filtered.txt;
    if [ ! -d $yacht_repo_loc/pathogen_detection_reference/virus_genomes ]; then
        mkdir $yacht_repo_loc/pathogen_detection_reference/virus_genomes;
        less $yacht_repo_loc/pathogen_detection_reference/virus_taxids_filtered.txt | cut -f 2 | parallel -j $cpu_num --link download_virus_assembly $yacht_repo_loc {};
    fi

    cd $yacht_repo_loc/pathogen_detection_reference;
    # download fungal genomes
    taxonkit list --ids 4751 --indent "" -r | grep 'species' | cut -d" " -f 1 | sort -u > fungi_taxids.txt
    taxonkit lineage -c -r -L fungi_taxids.txt > fungi_taxids_filtered.txt
    if [ ! -d $yacht_repo_loc/pathogen_detection_reference/fungi_genomes ]; then
        mkdir $yacht_repo_loc/pathogen_detection_reference/fungi_genomes;
        less $yacht_repo_loc/pathogen_detection_reference/fungi_taxids_filtered.txt | cut -f 2 | parallel -j $cpu_num --link download_fungi_assembly $yacht_repo_loc {};
    fi

    ## combine the microbe, virus and fungi genomes
    # find a representative genome for each species
    find $yacht_repo_loc/pathogen_detection_reference/microbe_genomes -name '*.fna' | sed 's/microbe_genomes\///' | sed 's/\.fna//' >> all_taxids.txt
    find $yacht_repo_loc/pathogen_detection_reference/virus_genomes -name '*.fna' | sed 's/virus_genomes\///' | sed 's/\.fna//' >> all_taxids.txt
    find $yacht_repo_loc/pathogen_detection_reference/fungi_genomes -name '*.fna' | sed 's/fungi_genomes\///' | sed 's/\.fna//' >> all_taxids.txt
    python $yacht_reproducibles_dir/benchmark/scripts/python_scripts/find_unique_spcies.py --taxids all_taxids.txt --target_taxids target_organisms.txt --output unique_species_taxids.tsv
    less all_taxids.txt | sort -u >> taxids.tmp
    less unique_species_taxids.tsv | cut -f 1 | sed '1d' | sort -u >> taxids.tmp
    less taxids.tmp | sort | uniq -u > delete_taxids.txt
    ## manually add some taxids (that are pretty close to target organisms)
    echo 2878546 >> delete_taxids.txt

    cd $yacht_repo_loc/pathogen_detection_reference;
    if [ ! -d $yacht_repo_loc/pathogen_detection_reference/all_genomes ]; then
        mkdir $yacht_repo_loc/pathogen_detection_reference/all_genomes;
        cd $yacht_repo_loc/pathogen_detection_reference/all_genomes;
        find $yacht_repo_loc/pathogen_detection_reference/microbe_genomes -name '*.fna' | xargs -I {} ln -s {} .;
        find $yacht_repo_loc/pathogen_detection_reference/virus_genomes -name '*.fna' | xargs -I {} ln -s {} .;
        find $yacht_repo_loc/pathogen_detection_reference/fungi_genomes -name '*.fna' | xargs -I {} ln -s {} .;
    fi
    # delete the genomes that are not representative
    less $yacht_repo_loc/pathogen_detection_reference/delete_taxids.txt | parallel -j $cpu_num rm $yacht_repo_loc/pathogen_detection_reference/all_genomes/{}.fna

    # build pathogen detection reference database for YACHT
    cd $yacht_repo_loc/pathogen_detection_reference;
    find ./all_genomes -name "*.fna" | awk -F"/" 'BEGIN{print "name,genome_filename,protein_filename"; OFS=","} {name=$3; sub(/\.fna$/, "", name); print name, $0 ","}' > datasets.csv
    sourmash sketch fromfile datasets.csv -p k=31,scaled=100,dna,abund -o customized_pathogen_detection_db_k31_scale100.zip
fi


## Crate the reference dictionary matrixes
# create a reference dictionary matrix for pathogen detection database
echo "Creating a reference dictionary matrix for pathogen detection database"
cd $yacht_repo_loc/pathogen_detection_reference;
if [ ! -f $yacht_repo_loc/pathogen_detection_reference/pathogen_detection_ani_thresh_0.995_config.json ]; then
    python $yacht_repo_loc/make_training_data_from_sketches.py --ref_file $yacht_repo_loc/pathogen_detection_reference/customized_pathogen_detection_db_k31_scale100.zip --ksize 31 --out_prefix 'pathogen_detection_ani_thresh_0.995' --ani_thresh 0.995
fi

## run YACHT on CAMI2 data
# create a "results" directory for the output
if [ ! -d $real_world_experiment_dir/results ]; then
    mkdir $real_world_experiment_dir/results
fi

# run YACHAT on the real-world data
echo "Running YACHT on the real-world data"
if [ ! -d $real_world_experiment_dir/results/YACHT_results ]; then
    mkdir $real_world_experiment_dir/results/YACHT_results;
    ## create sketches of samples
    fastq_list=`ls $real_world_experiment_dir/data/fastq`;
    if [ ! -d $real_world_experiment_dir/data/sketches ]; then
        mkdir $real_world_experiment_dir/data/sketches;
        cd $real_world_experiment_dir/data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=100,abund -o {}.sig.zip $real_world_experiment_dir/data/fastq/{}/*.fastq.gz ::: $fastq_list;
    fi
    ## run YACHT on the samples
    samples=`ls $real_world_experiment_dir/data/sketches | while read i;do filename=$(echo $i | sed 's/.sig.zip//'); echo $filename;done`
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --keep_raw --json $yacht_repo_loc/pathogen_detection_reference/pathogen_detection_ani_thresh_0.995_config.json --sample_file $real_world_experiment_dir/data/sketches/{}.fastq.gz.sig.zip --significance 0.99 --outdir $real_world_experiment_dir/results/YACHT_results --out_filename {}.xlsx ::: $samples;
fi
