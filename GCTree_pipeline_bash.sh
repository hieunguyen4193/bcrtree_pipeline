# install GCTree.
# see instruction to install GCTree https://matsen.group/gctree/install.html
# probably the best way is to use conda or pip install gctree
# install PHYLIP: conda install -c bioconda phylip

# example input args
# samplesheet="SampleSheet_FASTA.csv";
# deduplicate_src="./deduplicated.py"
# outputdir="./output/GCTree"

while getopts "i:o:d:" opt; do
  case ${opt} in
    i )
      samplesheet=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    d )
      deduplicate_src=$OPTARG
      ;;
    
    \? )
      echo "Usage: cmd [-i] samplesheet [-o] outputdir [-d] deduplicate_src"
      exit 1
      ;;
  esac
done

mkdir -p ${outputdir}
# if you don't use nextflow as a pipeline orchestration engine, run this bash script. 
while IFS=',' read -r sample_id fasta; do \
    mkdir -p ${outputdir}/${sample_id};
    echo -e "working on FASTA file " $sample_id " at " $fasta;
    fasta_name=$(echo ${fasta} | xargs -n 1 basename);

    cat ${fasta} | sed "s/>.*|Abundance:\\([0-9]\\+\\)/>\\1/" > ${outputdir}/${sample_id}/${fasta_name}.modified.fa
    python ${deduplicate_src} \
        --input ${outputdir}/${sample_id}/${fasta_name}.modified.fa \
        --root GL \
        --frame 0 \
        --id_abundances \
        --output_name ${sample_id} \
        --output ${outputdir}/${sample_id}

    mkconfig --quick ${outputdir}/${sample_id}/${sample_id}.phylip dnapars > ${outputdir}/${sample_id}/${sample_id}_dnapars.cfg
    dnapars < ${outputdir}/${sample_id}/${sample_id}_dnapars.cfg > ${outputdir}/${sample_id}/${sample_id}_dnapars.log

    # since outfile and outtree are generated at the current working dir, we need to manually move it to the output dir. 
    mv outfile ${outputdir}/${sample_id};
    mv outtree ${outputdir}/${sample_id};

    export QT_QPA_PLATFORM=offscreen
    export XDG_RUNTIME_DIR=/tmp/runtime-runner
    export MPLBACKEND=agg

    gctree infer \
        --verbose \
        --root GL \
        --idlabel ${outputdir}/${sample_id}/outfile \
        ${outputdir}/${sample_id}/${sample_id}.abundance.csv \
        --outbase ${sample_id};
        # --frame 1 \
        
    mv ${sample_id}* ${outputdir}/${sample_id};

done < <(tail -n +2 ${samplesheet});

echo -e "*******************************";
echo -e "***** FINISHED *****";
echo -e "*******************************";

# ***** DO NOT RUN *****
# test read in samplesheet file with sample id and path to file fasta. 
# samplesheet="SampleSheet_FASTA.csv";
# # if you don't use nextflow as a pipeline orchestration engine, run this bash script. 
# while IFS=',' read -r sampleid fasta; do \
#     echo $sampleid 
#     echo $fasta
# done < <(tail -n +2 ${samplesheet})
#####