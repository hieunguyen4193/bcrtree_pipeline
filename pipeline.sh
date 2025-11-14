# ***** run the mixcr pipeline to do alignment and BCR gene quantification.
# activate the conda env
# conda activate gctree

# export mixcr to current path
export PATH=/home/hieunguyen/CRC1382/src_2023/bcrtree_pipeline/mixcr:$PATH # change this to the path of your mixcr

fastqdir="/home/hieunguyen/CRC1382/src_2023/bcrtree_pipeline/input";
read_pattern="^N{35}(UMI:N{16})N{4:6}(R1:*)\^N{34}(R2:*)"; # to see more about tag pattern: https://mixcr.com/mixcr/reference/ref-tag-pattern/
outputdir="./output/mixcr_pipeline_output";
mkdir -p ${outputdir};

samplesheet="SampleSheet.csv";

while IFS=',' read -r sampleid fastq1 fastq2; do \
    if [ ! -f "${outputdir}/${sampleid}/${sampleid}.finished.txt" ]; then 
        echo -e "working on sample " ${sampleid}
        echo -e "path to FASTQ1 " ${fastq1}
        echo -e "path to FASTQ2 " ${fastq2}
        # run the mixcr pipeline
        bash mixcr_pipeline.sh -i ${sampleid} -o ${outputdir}  -f ${fastq1} -r ${fastq2} -p ${read_pattern}
        touch ${outputdir}/${sampleid}/${sampleid}.finished.txt
    else 
        echo -e "Sample " ${sampleid} " already processed, skipping..."
    fi
done < <(tail -n +2 ${samplesheet})

# ***** DO NOT RUN, THIS REQUIRES NEXTFLOW
# path_to_samplesheet=$1;
# output=$2;
# deduplicate_src=$3;
# path_to_work=$4;

# nextflow run GCTree_pipeline.nf \
# --samplesheet ${path_to_samplesheet} \
# --deduplicate_src ${deduplicate_src} \
# --output ${output} -resume -w ${path_to_work}

# ***** run "preprocessing" script to generate FASTA files. 
# each FASTA file contains all sequences used in constructing TREE. 
# we applied the criteria: clones must have at least 85% similarity in CDR3 sequences and same V-J genes 
# to be grouped into the same clone. 

# input=${outputdir}; # input to this = output of mixcr
# fasta="./output/fasta_for_GCTree_pipeline";
# processed_output="./output/processed_output"
# PROJECT="GCTree_mixcr_test";
# re_define_clone_cluster=TRUE
# rerun=FALSE
# savefile=TRUE
# verbose=TRUE
# save_fasta=TRUE
# define_clone_clusters=FALSE
# ref_gene="IMGT"
# mouse_id="testMouse"

# Rscript prepare_input_for_GCTree.R \
#     --input ${input} \
#     --fasta ${fasta} \
#     --processed_output ${processed_output} \
#     --PROJECT ${PROJECT} \
#     --re_define_clone_cluster ${re_define_clone_cluster} \
#     --rerun ${rerun} \
#     --savefile ${savefile} \
#     --verbose ${verbose} \
#     --save_fasta ${save_fasta} \
#     --define_clone_cluster ${define_clone_clusters} \
#     --ref_gene ${ref_gene} \
#     --mouse_id ${mouse_id};
# the above rscript also generates a file "SampleSheet_FASTA.csv"

# ***** GCTree -> generate tree
# the file "SampleSheet_FASTA.csv" contains filename and path of all FASTA files. 
# 1 fasta file --> 1 tree (if possible, if a parsimony tree exists)
bash GCTree_pipeline_bash.sh -i "SampleSheet_FASTA.csv" -o "./output/GCTree" -d "./deduplicated.py"