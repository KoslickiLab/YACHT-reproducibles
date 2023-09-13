#! /bin/bash
# run benchmarking experiments based on CAMI2 data for YACHT algorithm
# Usage:    

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run_YACHT.sh <yacht_repo_loc> <benchmark_dir> <cpu_num>"
    exit 1
fi
yacht_repo_loc=$1
benchmark_dir=$2
cpu_num=$3

## Run YACHT on CAMI2 data (e.g., Rhizosphere challenge, Clinical pathogen detection challenge, Challenge Marine Dataset, Strain Madness Dataset)

# Build a NCBI taxonomy database used for CAMI2 Challenge
cd $yacht_repo_loc
if [ ! -d $yacht_repo_loc/ncbi_reference ]; then
    mkdir $yacht_repo_loc/ncbi_reference;
    wget -p $yacht_repo_loc/ncbi_reference https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_2_DATABASES/ncbi_taxonomy.tar;

    # copy names.dmp, nodes.dmp, delnodes.dmp and merged.dmp to data directory: $HOME/.taxonkit
    cd $HOME/.taxonkit;
    ln -s $yacht_repo_loc/ncbi_reference/ncbi_taxonomy.tar;
    tar xvf ncbi_taxonomy.tar;
    mv ./ncbi_taxonomy/taxdump.tar.gz ./;
    tar zxvf taxdump.tar.gz;
    wget -p $yacht_repo_loc/ncbi_reference https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_2_DATABASES/RefSeq_genomic_20190108.tar;
    cd $yacht_repo_loc/ncbi_reference;
    mkdir RefSeq_genomic_20190108;
    tar xvf RefSeq_genomic_20190108.tar -C RefSeq_genomic_20190108;

    # find the mapping between RefSeq and taxid
    find ./RefSeq_genomic_20190108/ -type f | parallel -j $cpu_num "refseq=\$(echo {} | sed 's/\.\/RefSeq_genomic_20190108\///' | cut -d'_' -f 1-2); taxid=\$(datasets summary genome accession \$refseq | jq '.reports[0].organism.tax_id'); echo \$refseq,\$taxid" > refseq_taxid_mapping.csv;

    # build database metadata file
    python $benchmark_dir/scripts/python_scripts/build_ncbi_database_metadata.py --metadata $benchmark_dir/CAMI_data/rhizosphere_data/genomes/combined_metadata.tsv $benchmark_dir/CAMI_data/marine_data/genomes/combined_metadata.tsv $benchmark_dir/CAMI_data/strain_madness_data/genomes/combined_metadata.tsv --refseq_taxid_mapping $yacht_repo_loc/ncbi_reference/refseq_taxid_mapping.csv;

    # build NCBI reference database for YACHT
    less ncbi_database_metadata.tsv | sed '1d' | awk -F"\t" 'BEGIN{print "name,genome_filename,protein_filename"; OFS=","}{print $1,$4,""}' > datasets.csv
    sourmash sketch fromfile datasets.csv -p k=31,scaled=100,dna,abund -o customized_ncbi_db_k31_scale100.zip
fi

# Build a pathogen database (data collected from BVBRC) used for CAMI2 Challenge
cd $yacht_repo_loc
if [ ! -d $yacht_repo_loc/pathogen_detection_reference ]; then
    mkdir $yacht_repo_loc/pathogen_detection_reference;
    cd $yacht_repo_loc/pathogen_detection_reference;
    # download non-viral genomes from BVBRC
    wget -O $yacht_repo_loc/pathogen_detection_reference/microbe_genome_metadata.tsv ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata;
    ## add quotation mark to avoid float matching issue
    # less microbe_genome_metadata.tsv  | awk 'BEGIN {FS=OFS="\t"} NR>1 {$1="'"'"'"$1"'"'"'"} 1' > a
    # mv a microbe_genome_metadata.tsv
    awk 'BEGIN{IGNORECASE=1;FS="\t";OFS="\t"} $46 ~ /human; |Homo sapiens/ {print $1,$2,$4,$5,$18,$19,$20,$46,$64}' microbe_genome_metadata.tsv | cut -f 3 | sort -u > microbe_taxids.txt
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
    fi

    cd $yacht_repo_loc/pathogen_detection_reference;
    # download viral genomes from BVBRC
    wget -O $yacht_repo_loc/pathogen_detection_reference/virus_genome_metadata.tsv ftp://ftp.bvbrc.org/viruses/genome_metadata;
    # less virus_genome_metadata.tsv  | awk 'BEGIN {FS=OFS="\t"} NR>1 {$1="'"'"'"$1"'"'"'"} 1' > a
    # mv a virus_genome_metadata.tsv
    awk 'BEGIN{IGNORECASE=1;FS="\t";OFS="\t"} $46 ~ /human; |Homo sapiens/ {print $1,$2,$4,$5,$18,$19,$20,$46,$64}' virus_genome_metadata.tsv | cut -f 3 | sort -u > virus_taxids.txt
    taxonkit lineage -c -r -L virus_taxids.txt | awk -F'\t' '$3 ~ /rank|species|strain/ {print $0}' > virus_taxids_filtered.txt;
    if [ ! -d $yacht_repo_loc/pathogen_detection_reference/virus_genomes ]; then
        mkdir $yacht_repo_loc/pathogen_detection_reference/virus_genomes;
        cd $yacht_repo_loc/pathogen_detection_reference/virus_genomes;
        less $yacht_repo_loc/pathogen_detection_reference/virus_taxids_filtered.txt | cut -f 2 | xargs -I {} -P $cpu_num sh -c 'mkdir $0/pathogen_detection_reference/virus_genomes/{}; cd $0/pathogen_detection_reference/virus_genomes/{}; ncbi-genome-download -s refseq -F "fasta" -l "complete" -t {} viral; ncbi-genome-download -s genbank -F "fasta" -l "complete" -t {} viral' $yacht_repo_loc;
        cd $yacht_repo_loc/pathogen_detection_reference/virus_genomes;
        for taxid in *; do
            filename=`find ./$taxid -type f -name "*.fna.gz" | sort -t'/' -k4 -r | head -1`
            if [ -n "$filename" ]; then
                echo "$taxid" >> $yacht_repo_loc/pathogen_detection_reference/virus_check_list.txt;
                mv $filename ${taxid}.fna.gz;
                gunzip ${taxid}.fna.gz;    
            fi
            rm -rf $taxid;
        done
    fi

    ## combine the microbe and virus genomes
    cd $yacht_repo_loc/pathogen_detection_reference;
    if [ ! -d $yacht_repo_loc/pathogen_detection_reference/all_genomes ]; then
        mkdir $yacht_repo_loc/pathogen_detection_reference/all_genomes;
        cd $yacht_repo_loc/pathogen_detection_reference/all_genomes;
        ln -s $yacht_repo_loc/pathogen_detection_reference/microbe_genomes/* .;
        ln -s $yacht_repo_loc/pathogen_detection_reference/virus_genomes/* .;
    fi

    # build pathogen detection reference database for YACHT
    cd $yacht_repo_loc/pathogen_detection_reference;
    find ./all_genomes -name "*.fna" | awk -F"/" 'BEGIN{print "name,genome_filename,protein_filename"; OFS=","} {name=$3; sub(/\.fna$/, "", name); print name, $0 ","}' > datasets.csv
    sourmash sketch fromfile datasets.csv -p k=31,scaled=100,dna,abund -o customized_pathogen_detection_db_k31_scale100.zip
fi

## Crate the reference dictionary matrixes
# create a reference dictionary matrix for NCBI database
cd $yacht_repo_loc/ncbi_reference;
if [ ! -f $yacht_repo_loc/ncbi_reference/ncbi_ani_thresh_0.995_ref_matrix_processed.npz ]; then
    python $yacht_repo_loc/make_training_data_from_sketches.py --ref_file $yacht_repo_loc/ncbi_reference/customized_ncbi_db_k31_scale100.zip --ksize 31 --out_prefix 'ncbi_ani_thresh_0.995' --ani_thresh 0.995
fi


# create a reference dictionary matrix for pathogen detection database
cd $yacht_repo_loc/pathogen_detection_reference;
if [ ! -f $yacht_repo_loc/pathogen_detection_reference/pathogen_detection_ani_thresh_0.95_ref_matrix_processed.npz ]; then
    python $yacht_repo_loc/make_training_data_from_sketches.py --ref_file $yacht_repo_loc/pathogen_detection_reference/customized_pathogen_detection_db_k31_scale100.zip --ksize 31 --out_prefix 'pathogen_detection_ani_thresh_0.95' --ani_thresh 0.95
fi

## run YACHT on CAMI2 data
# create a "results" directory for the output
if [ ! -d $benchmark_dir/results ]; then
    mkdir $benchmark_dir/results
fi
if [ ! -d $benchmark_dir/results/YACHT_results ]; then
    mkdir $benchmark_dir/results/YACHT_results
fi

# run YACHAT on Rhizosphere challenge data
if [ ! -d $benchmark_dir/results/YACHT_results/rhizosphere_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/rhizosphere_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/rhizosphere_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/rhizosphere_data/sketches;
        cd $benchmark_dir/CAMI_data/rhizosphere_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=100,abund -o rhimgCAMI2_sample_{}.sig.zip $benchmark_dir/CAMI_data/rhizosphere_data/rhimgCAMI2_sample_{}.fq.gz ::: {0..20};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/ncbi_reference/ncbi_ani_thresh_0.995_ref_matrix_processed.npz --sample_file $benchmark_dir/CAMI_data/rhizosphere_data/sketches/rhimgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.995 --significance 0.99 --min_coverage 1 --outfile $benchmark_dir/results/YACHT_results/rhizosphere_data/rhimgCAMI2_sample_{}.csv ::: {0..20};
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format_ncbi.py --yacht_res $benchmark_dir/results/YACHT_results/rhizosphere_data/rhimgCAMI2_sample_{2}.csv --metadata_dir $yacht_repo_loc/ncbi_reference --min_coverage {1} --outfile $benchmark_dir/results/YACHT_results/rhizosphere_data/rhimgCAMI2_sample_{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: {0..20};
fi

# run YACHT on Clinical pathogen detection challenge data
if [ ! -d $benchmark_dir/results/YACHT_results/pathogen_detection_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/pathogen_detection_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/pathogen_detection_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/pathogen_detection_data/sketches;
        cd $benchmark_dir/CAMI_data/pathogen_detection_data/sketches;
        sourmash sketch dna -f -p k=31,scaled=100,abund -o patmgCAMI2.sig.zip $benchmark_dir/CAMI_data/pathogen_detection_data/*.fastq.gz;
    fi
    ## run YACHT on the samples
    python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/pathogen_detection_reference/pathogen_detection_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $benchmark_dir/CAMI_data/pathogen_detection_data/sketches/patmgCAMI2.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2.csv;
    ## convert the YACHT results to CAMI format
    # python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format.py --yacht_res $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2.csv --metadata_dir $yacht_repo_loc/ncbi_reference --outfile $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2_cami_format.profile;
fi

# run YACHT on Challenge Marine Dataset
if [ ! -d $benchmark_dir/results/YACHT_results/marine_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/marine_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/marine_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/marine_data/sketches;
        cd $benchmark_dir/CAMI_data/marine_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=100,abund -o marmgCAMI2_sample_{}.sig.zip $benchmark_dir/CAMI_data/marine_data/marmgCAMI2_sample_{}.fq.gz ::: {0..9};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/ncbi_reference/ncbi_ani_thresh_0.995_ref_matrix_processed.npz --sample_file $benchmark_dir/CAMI_data/marine_data/sketches/marmgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.995 --significance 0.99 --min_coverage 1 --outfile $benchmark_dir/results/YACHT_results/marine_data/marmgCAMI2_sample_{}.csv ::: {0..9};
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format_ncbi.py --yacht_res $benchmark_dir/results/YACHT_results/marine_data/marmgCAMI2_sample_{2}.csv --metadata_dir $yacht_repo_loc/ncbi_reference --min_coverage {1} --outfile $benchmark_dir/results/YACHT_results/marine_data/marmgCAMI2_sample_{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: {0..9};
fi

# run YACHT on Strain Madness Dataset
if [ ! -d $benchmark_dir/results/YACHT_results/strain_madness_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/strain_madness_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/strain_madness_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/strain_madness_data/sketches;
        cd $benchmark_dir/CAMI_data/strain_madness_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=100,abund -o strmgCAMI2_sample_{}.sig.zip $benchmark_dir/CAMI_data/strain_madness_data/strmgCAMI2_sample_{}.fq.gz ::: {0..99};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/ncbi_reference/ncbi_ani_thresh_0.995_ref_matrix_processed.npz --sample_file $benchmark_dir/CAMI_data/strain_madness_data/sketches/strmgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.995 --significance 0.99 --min_coverage 1 --outfile $benchmark_dir/results/YACHT_results/strain_madness_data/strmgCAMI2_sample_{}.csv ::: {0..99};
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format_ncbi.py --yacht_res $benchmark_dir/results/YACHT_results/strain_madness_data/strmgCAMI2_sample_{2}.csv --metadata_dir $yacht_repo_loc/ncbi_reference --min_coverage {1} --outfile $benchmark_dir/results/YACHT_results/strain_madness_data/strmgCAMI2_sample_{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: {0..99};
fi