while getopts "i:o:f:r:p:" opt; do
  case ${opt} in
    i )
      sampleid=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    f )
      fastq1=$OPTARG
      ;;
    r )
      fastq2=$OPTARG
      ;;
    p )
      read_pattern=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] sampleid [-o] outputdir [-f] fastq1 [-r] fastq2 [-p]  read_pattern"
      exit 1
      ;;
  esac
done

output=${outputdir}/${sampleid};
mkdir -p ${output};

# example read pattern, applied to Britta BCR data. 
# read_pattern="^N{35}(UMI:N{16})N{4:6}(R1:*)\^N{34}(R2:*)"

#####--------------------------------------------------#####
##### UPSTREAM PIPELINE
#####--------------------------------------------------#####

#### ALIGNMENT
echo -e "ALIGNMENT"
echo -e "using read pattern " ${read_pattern}
mixcr align \
-p generic-amplicon-with-umi \
--rna \
--library imgt.202312-3 \
--species mmu \
--report ${output}/${sampleid}.report.txt \
--json-report ${output}/${sampleid}.align.report.json \
--tag-pattern ${read_pattern} \
--floating-left-alignment-boundary \
--floating-right-alignment-boundary C \
$fastq1 \
$fastq2 \
${output}/${sampleid}.vdjca

##### Refine tags and sort
echo -e "Refine tags and sort"

mixcr refineTagsAndSort \
--report ${output}/${sampleid}.refine.report.txt \
--json-report ${output}/${sampleid}.refine.report.json \
${output}/${sampleid}.vdjca \
${output}/${sampleid}.refined.vdjca

##### Assemble
echo -e "Assemble"

mixcr assemble \
--report ${output}/${sampleid}.assemble.report.txt \
--json-report ${output}/${sampleid}.assemble.report.json \
-OassemblingFeatures={FR1Begin:FR4End} \
-OseparateByC=true \
${output}/${sampleid}.refined.vdjca \
${output}/${sampleid}.clns

##### Export CLONES
echo -e "EXRPOT CLONES"
mixcr exportClones -uniqueTagCount Molecule -count \
${output}/${sampleid}.clns ${output}/${sampleid}.tsv

##### Export QUALITY CONTROLS
echo -e "Export QC"
mixcr exportQc align ${output}/*.clns ${output}/alignQc.pdf

#####--------------------------------------------------#####
##### TREES
#####--------------------------------------------------#####
mixcr findAlleles \
--report ${output}/${sampleid}.findAlleles.report.txt \
--json-report ${output}/${sampleid}.findAlleles.report.json \
--export-alleles-mutations ${output}/${sampleid}_alleles.tsv \
--export-library ${output}/${sampleid}_alleles.json \
--output-template {file_dir_path}/{file_name}.reassigned.clns \
${output}/${sampleid}.clns

mixcr exportClones -uniqueTagCount Molecule -count ${output}/${sampleid}.reassigned.clns ${output}/${sampleid}.reassigned.tsv

# EOF
# to install MIXCR: https://mixcr.com/mixcr/getting-started/installation/
# download mixcr 4.6: https://github.com/milaboratory/mixcr/releases/download/v4.6.0/mixcr-4.6.0.zip
# tested with version 4.3.2
# to get the license: https://platforma.bio/getlicense?_gl=1*1dx2w72*_gcl_au*MjE0MzIwOTAyLjE3NjI4NzQwNjM.*_ga*MTMxODQwNTg5NC4xNzYyODc0MDYz*_ga_7K1RPNHREB*czE3NjI4NzQwNjMkbzEkZzEkdDE3NjI4NzQwNzIkajUxJGwwJGgw&product=mixcr
# mixcr activate-license

# to dowload imgt  reference genome: https://mixcr.com/mixcr/guides/external-libraries/
# download from imgt github: https://github.com/repseqio/library-imgt/releases