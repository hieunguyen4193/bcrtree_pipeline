fastqdir="/home/hieunguyen/CRC1382/src_2023/bcrtree_pipeline/input"

# path_to_samplesheet=$1;
# output=$2;
# deduplicate_src=$3;
# path_to_work=$4;

nextflow run GCtree_pipeline_input_SampleSheet.nf \
--samplesheet ${path_to_samplesheet} \
--deduplicate_src ${deduplicate_src} \
--output ${output} -resume -w ${path_to_work}