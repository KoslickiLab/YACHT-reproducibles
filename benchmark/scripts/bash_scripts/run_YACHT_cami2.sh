#! /bin/bash
# run benchmarking experiments based on CAMI2 data for YACHT algorithm  

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run_YACHT_cami2.sh <yacht_reproducibles_dir> <cpu_num>"
    exit 1
fi

yacht_reproducibles_dir=$1
yacht_repo_loc=${yacht_reproducibles_dir}/YACHT
benchmark_dir=${yacht_reproducibles_dir}/benchmark
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



## Run YACHT on CAMI2 data (e.g., Rhizosphere challenge, Clinical pathogen detection challenge, Challenge Marine Dataset, Strain Madness Dataset)
# Build a NCBI taxonomy database used for CAMI2 Challenge
echo "Building a NCBI taxonomy database used for CAMI2 Challenge"
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
    echo 2878546 2785247 >> delete_taxids.txt

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

    ## group the genomes with the same MD5 hash (they are the same genome so we only keep one copy)
    python $yacht_reproducibles_dir/benchmark/scripts/python_scripts/group_same_genomes.py --sourmash_db customized_pathogen_detection_db_k31_scale100.zip;
fi

## Crate the reference dictionary matrixes
# create a reference dictionary matrix for NCBI database
echo "Creating a reference dictionary matrix for NCBI database"
cd $yacht_repo_loc/ncbi_reference;
if [ ! -f $yacht_repo_loc/ncbi_reference/ncbi_ani_thresh_0.995_config.json ]; then
    python $yacht_repo_loc/make_training_data_from_sketches.py --ref_file $yacht_repo_loc/ncbi_reference/customized_ncbi_db_k31_scale100.zip --ksize 31 --num_threads 200 --ani_thresh 0.995 --prefix 'ncbi_ani_thresh_0.995' --outdir ./ &
fi


# create a reference dictionary matrix for pathogen detection database
echo "Creating a reference dictionary matrix for pathogen detection database"
cd $yacht_repo_loc/pathogen_detection_reference;
if [ ! -f $yacht_repo_loc/pathogen_detection_reference/pathogen_detection_ani_thresh_0.995_config.json ]; then
    python $yacht_repo_loc/make_training_data_from_sketches.py --ref_file $yacht_repo_loc/pathogen_detection_reference/customized_pathogen_detection_db_k31_scale100_filtered.zip --ksize 31 --num_threads 200 --ani_thresh 0.995 --prefix 'pathogen_detection_ani_thresh_0.995' --outdir ./
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
echo "Running YACHT on Rhizosphere challenge data"
if [ ! -d $benchmark_dir/results/YACHT_results/rhizosphere_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/rhizosphere_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/rhizosphere_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/rhizosphere_data/sketches;
    fi
    cd $benchmark_dir/CAMI_data/rhizosphere_data/sketches;
    samples=`ls $benchmark_dir/CAMI_data/rhizosphere_data/*.fq.gz | awk -F"/" '{print $NF}' | sed 's/\.fq\.gz//'`;
    for sample_name in $samples; do
        if [ ! -d $benchmark_dir/CAMI_data/rhizosphere_data/sketches/${sample_name} ]; then
            mkdir $benchmark_dir/CAMI_data/rhizosphere_data/sketches/${sample_name};
        fi
    done
    parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=100,abund -o $benchmark_dir/CAMI_data/rhizosphere_data/sketches/{}/{}.sig.zip $benchmark_dir/CAMI_data/rhizosphere_data/{}.fq.gz ::: $samples;
    # run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --keep_raw --json $yacht_repo_loc/ncbi_reference/ncbi_ani_thresh_0.995_config.json --sample_file $benchmark_dir/CAMI_data/rhizosphere_data/sketches/{}/{}.sig.zip --significance 0.99 --num_threads 100 --out $benchmark_dir/results/YACHT_results/rhizosphere_data/{}.xlsx ::: $samples;
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format_ncbi.py --yacht_res $benchmark_dir/results/YACHT_results/rhizosphere_data/{2}.xlsx --metadata_dir $yacht_repo_loc/ncbi_reference --min_coverage {1} --outfile $benchmark_dir/results/YACHT_results/rhizosphere_data/{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: $samples;
fi


# run YACHT on Clinical pathogen detection challenge data
echo "Running YACHT on Clinical pathogen detection challenge data"
if [ ! -d $benchmark_dir/results/YACHT_results/pathogen_detection_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/pathogen_detection_data;
    echo "Unzipping the input file..."
    if [ ! -e "${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq" ]; then
        gunzip -c ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq.gz > ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq
        metagenome_1=${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq
    else
        metagenome_1=${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R1.fastq
    fi
    if [ ! -e "${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq" ]; then
        gunzip -c ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq.gz > ${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq
        metagenome_2=${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq
    else
        metagenome_2=${benchmark_dir}/CAMI_data/pathogen_detection_data/patmg_CAMI2_short_read_R2.fastq
    fi
    # combine two fastq files into one
    echo "Combining two fastq files into one..."
    if [ ! -e "${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq" ]; then
        cat ${metagenome_1} ${metagenome_2} > ${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq
        metaphlan_in=${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq
    else
        metaphlan_in=${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq
    fi
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/pathogen_detection_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/pathogen_detection_data/sketches;

    fi
    cd $benchmark_dir/CAMI_data/pathogen_detection_data/sketches;
    sourmash sketch dna -f -p k=31,scaled=100,abund -o patmgCAMI2.sig.zip $metaphlan_in;
    ## run YACHT on the samples
    python $yacht_repo_loc/run_YACHT.py --keep_raw --json $yacht_repo_loc/pathogen_detection_reference/pathogen_detection_ani_thresh_0.995_config.json --sample_file $benchmark_dir/CAMI_data/pathogen_detection_data/sketches/patmgCAMI2.sig.zip --num_threads 100 --significance 0.99  --out $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2_0.995.xlsx;
    ## convert the YACHT results to CAMI format
    # python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format.py --yacht_res $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2.csv --metadata_dir $yacht_repo_loc/ncbi_reference --outfile $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2_cami_format.profile;
fi
metaphlan_in=${benchmark_dir}/CAMI_data/pathogen_detection_data/combined.fastq
cd $benchmark_dir/CAMI_data/pathogen_detection_data/sketches;
sourmash sketch dna -f -p k=31,scaled=100,abund -o patmgCAMI2.sig.zip $metaphlan_in;
python $yacht_repo_loc/run_YACHT.py --keep_raw --json $yacht_repo_loc/pathogen_detection_reference/pathogen_detection_ani_thresh_0.995_config.json --sample_file $benchmark_dir/CAMI_data/pathogen_detection_data/sketches/patmgCAMI2.sig.zip --num_threads 100 --significance 0.99  --out $benchmark_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2_0.995.xlsx;


# run YACHT on Challenge Marine Dataset
echo "Running YACHT on Challenge Marine Dataset"
if [ ! -d $benchmark_dir/results/YACHT_results/marine_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/marine_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/marine_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/marine_data/sketches;
    fi
    cd $benchmark_dir/CAMI_data/marine_data/sketches;
    samples=`ls $benchmark_dir/CAMI_data/marine_data/*.fq.gz | awk -F"/" '{print $NF}' | sed 's/\.fq\.gz//'`;
    for sample_name in $samples; do
        if [ ! -d $benchmark_dir/CAMI_data/marine_data/sketches/${sample_name} ]; then
            mkdir $benchmark_dir/CAMI_data/marine_data/sketches/${sample_name};
        fi
    done
    parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=100,abund -o $benchmark_dir/CAMI_data/rhizosphere_data/{}/{}.sig.zip $benchmark_dir/CAMI_data/marine_data/{}.fq.gz ::: $samples;
    # run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --keep_raw --json $yacht_repo_loc/ncbi_reference/ncbi_ani_thresh_0.995_config.json --sample_file $benchmark_dir/CAMI_data/marine_data/sketches/{}/{}.sig.zip --significance 0.99 --num_threads 100 --out $benchmark_dir/results/YACHT_results/marine_data/{}.xlsx ::: $samples;
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format_ncbi.py --yacht_res $benchmark_dir/results/YACHT_results/marine_data/{2}.xlsx --metadata_dir $yacht_repo_loc/ncbi_reference --min_coverage {1} --outfile $benchmark_dir/results/YACHT_results/marine_data/{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: $samples;
fi

# run YACHT on Strain Madness Dataset
echo "Running YACHT on Strain Madness Dataset"
if [ ! -d $benchmark_dir/results/YACHT_results/strain_madness_data ]; then
    mkdir $benchmark_dir/results/YACHT_results/strain_madness_data;
    ## create sketches of samples
    if [ ! -d $benchmark_dir/CAMI_data/strain_madness_data/sketches ]; then
        mkdir $benchmark_dir/CAMI_data/strain_madness_data/sketches;
    fi
    cd $benchmark_dir/CAMI_data/strain_madness_data/sketches;
    samples=`ls $benchmark_dir/CAMI_data/strain_madness_data/*.fq.gz | awk -F"/" '{print $NF}' | sed 's/\.fq\.gz//'`;
    for sample_name in $samples; do
        if [ ! -d $benchmark_dir/CAMI_data/strain_madness_data/sketches/${sample_name} ]; then
            mkdir $benchmark_dir/CAMI_data/strain_madness_data/sketches/${sample_name};
        fi
    done
    parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=100,abund -o $benchmark_dir/CAMI_data/strain_madness_data/{}/{}.sig.zip $benchmark_dir/CAMI_data/strain_madness_data/{}.fq.gz ::: $samples;
    # run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --keep_raw --json $yacht_repo_loc/ncbi_reference/ncbi_ani_thresh_0.995_config.json --sample_file $benchmark_dir/CAMI_data/strain_madness_data/sketches/{}/{}.sig.zip --significance 0.99 --num_threads 100 --out $benchmark_dir/results/YACHT_results/strain_madness_data/{}.xlsx ::: $samples;
    ## convert the YACHT results to CAMI format
    parallel -j $cpu_num python $benchmark_dir/scripts/python_scripts/convert_to_CAMI_format_ncbi.py --yacht_res $benchmark_dir/results/YACHT_results/strain_madness_data/{2}.xlsx --metadata_dir $yacht_repo_loc/ncbi_reference --min_coverage {1} --outfile $benchmark_dir/results/YACHT_results/strain_madness_data/{2}_cami_format_coverage{1}.profile ::: 1 0.5 0.1 0.05 0.01 ::: $samples;
fi