# ***** run the mixcr pipeline to do alignment and BCR gene quantification.
export PATH=/home/hieunguyen/CRC1382/src_2023/bcrtree_pipeline/mixcr:$PATH
fastqdir="/home/hieunguyen/CRC1382/src_2023/bcrtree_pipeline/input";
read_pattern="^N{35}(UMI:N{16})N{4:6}(R1:*)\^N{34}(R2:*)"; # to see more about tag pattern: https://mixcr.com/mixcr/reference/ref-tag-pattern/
outputdir="./output/mixcr_pipeline_output";
mkdir -p ${outputdir};

samplesheet="SampleSheet.csv";

while IFS=',' read -r sampleid fastq1 fastq2; do
    echo -e "working on sample " ${sampleid}
    echo -e "path to FASTQ1 " ${fastq1}
    echo -e "path to FASTQ2 " ${fastq2}
    # run the mixcr pipeline
    bash mixcr_pipeline.sh -i ${sampleid} -o ${outputdir}  -f ${fastq1} -r ${fastq2} -p ${read_pattern}
done < <(tail -n +2 ${samplesheet})

# *****
# path_to_samplesheet=$1;
# output=$2;
# deduplicate_src=$3;
# path_to_work=$4;

# nextflow run GCTree_pipeline.nf \
# --samplesheet ${path_to_samplesheet} \
# --deduplicate_src ${deduplicate_src} \
# --output ${output} -resume -w ${path_to_work}