# General information 
This pipeline takes `FASTQ` files and do B-cell repertoire sequence alignment and V-J gene quantification. After that, we generate BCR evolutional trees using GCTree algorithm, see the paper [Using Genotype Abundance to Improve Phylogenetic Inference](https://academic.oup.com/mbe/article/35/5/1253/4893244).

# Preparation
We need the following tools/packages to run the pipeline: 
## Install `MiXCR`
To run this pipeline, you need to install: 
- Download `MiXCR` package from their github repo: 
```sh
wget https://github.com/milaboratory/mixcr/releases/download/v4.6.0/mixcr-4.6.0.zip
```
- Unzip the file and follow the instruction in [installtion guide](https://mixcr.com/mixcr/getting-started/installation/)
- Export `MiXCR` to your current path (to make it executable):
```sh
export PATH=/where/you/store/mixcr:$PATH
```

## Install `GCTree`
Install the package `GCTree` in a `conda` environment:
```sh 
# Run this script one line at a time
conda create -n gctree python=3.9 -y
conda activate gctree
pip install gctree
conda install -c bioconda phylip -y
```

## R
`R` is also required. I've already added script to automatically check and install missing packages in the R environment. In case there is still some packages missing, please install them manually.

```R 
# a code snippet to install R packages, taken from "prepare_input_for_GCTree.R"
new.pkgs <- c("stringr", "tidyverse", "dplyr", 
              "ggplot2", "stringdist", "igraph",
              "this.path", "argparse")
new.bioc.pkgs <- c("Biostrings", "msa")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}
for (pkg in new.bioc.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    BiocManager::install(pkg)
  }
}
```

# Run the pipeline
The pipeline is written in the script `pipeline.sh`. I'm going to break down the script and explain its function:

We need to export the `MiXCR` tool to make it executable:
```sh
export PATH=/home/hieunguyen/CRC1382/src_2023/bcrtree_pipeline/mixcr:$PATH # change this to the path of your mixcr
```

```sh
fastqdir="/home/hieunguyen/CRC1382/src_2023/bcrtree_pipeline/input";
read_pattern="^N{35}(UMI:N{16})N{4:6}(R1:*)\^N{34}(R2:*)"; # to see more about tag pattern: https://mixcr.com/mixcr/reference/ref-tag-pattern/
outputdir="./output/mixcr_pipeline_output";
mkdir -p ${outputdir};

samplesheet="SampleSheet.csv";
```

The variable `fastqdir` defines path to the directory containing all input `FASTQ` files. For example, I store all my `FASTQ` files in the folder 
```sh
"/home/hieunguyen/CRC1382/src_2023/bcrtree_pipeline/input"
```
You can download these `FASTQ` files to replicate my analysis from the One Drive link. 

`read_pattern` defines the structure of the read - the library. For example, here we have reads which have 35bp of barcode, a UMI of 16bp long, a random sequence barcode of 4-6bp long, then read R1. This is followed by a 34bp barcode again and then the rest is read R2. You need to modify this `read_pattern` variable to match your library. 

All output are stored in the path specified at `outputdir`. In this example, I export all outputs to a folder named `output` in the same working directory. 

The file `SampleSheet.csv` is where we specify our input `FASTQ` files. This is a `.csv` file with 3 columns: `SampleID`, `fastq1` and `fastq2` (pair-end sequencing data). Name of the samples are written in the `SampleID` column and they will be used in all downstream tasks. The column `fastq1` and `fastq2` specify paths to the `R1.fastq` and `R2.fastq`, respectively. 

We start running the `MiXCR` pipeline: 

```sh
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
```

This script loops through each line of the `SampleSheet.csv` file, read in the `SampleID`, the `FASTQ` files' paths `fastq1, fastq2`. We execute the `MiXCR` pipeline via the command 
```sh
bash mixcr_pipeline.sh -i ${sampleid} -o ${outputdir}  -f ${fastq1} -r ${fastq2} -p ${read_pattern}
```

When finished, we generate an empty file `touch ${outputdir}/${sampleid}/${sampleid}.finished.txt` to mark that the process is finished. If you wanted to re-run the `MiXCR` for the same sample, you need to delete this file. Otherwise the pipeline will not be execute due to the `if-else` structure `if [ ! -f "${outputdir}/${sampleid}/${sampleid}.finished.txt" ]; then ...`.

Now we move to the next part: prepare the `MiXCR` output for `GCTree` tree inference. 
```sh
input=${outputdir}; # input to this = output of mixcr
fasta="./output/fasta_for_GCTree_pipeline";
processed_output="./output/processed_output"
PROJECT="GCTree_mixcr_test";
re_define_clone_cluster=TRUE
rerun=FALSE
savefile=TRUE
verbose=TRUE
save_fasta=TRUE
define_clone_clusters=FALSE
ref_gene="IMGT"
mouse_id="testMouse"

Rscript prepare_input_for_GCTree.R \
    --input ${input} \
    --fasta ${fasta} \
    --processed_output ${processed_output} \
    --PROJECT ${PROJECT} \
    --re_define_clone_cluster ${re_define_clone_cluster} \
    --rerun ${rerun} \
    --savefile ${savefile} \
    --verbose ${verbose} \
    --save_fasta ${save_fasta} \
    --define_clone_cluster ${define_clone_clusters} \
    --ref_gene ${ref_gene} \
    --mouse_id ${mouse_id};
```
- We use the output of `MiXCR` as input to the R-script `prepare_input_for_GCTree.R`: `input=${outputdir}`.
- This R-script will generate `FASTA` files. Each line in the `FASTA` file is a BCR sequence. The name of the sequence contains information on V-J genes, amino acid sequence, length of sequence, sample ID. Most importantly, it contains the **abundance** of the sequence. The **abundance** is also used in `GCTree` inference algorithm. This script saves `FASTA` files to the path specified in variable `fasta`. Intermediate files are saved to path in the variable `processed_output`. 
- `PROJECT`: name of the project. 
- `mouseID`: name of the mouse. We collect all BCR sequences from all samples of the same mouse and generate the tree. All samples (MID) in the *example input (in One Drive)* are from the same mouse. 
- `ref_gene="IMGT"`: we use the IMGT reference genome database for BCR sequences. 
- For other parameters, just keep as default. 

Finally, we run the `GCTree` script to generate trees for each `FASTA` file generated in previous step

```sh 
bash GCTree_pipeline_bash.sh -i "SampleSheet_FASTA.csv" -o "./output/GCTree" -d "./deduplicated.py"
```

If there is no tree found from the input BCR sequences, the following error will be shown
```sh
/home/hieunguyen/miniconda3/envs/gctree/lib/python3.9/site-packages/gctree/phylip_parse.py:189: UserWarning: Parsed 0 from ./output/GCTree/IGHV1-47-01_IGHJ4-01_42_1/outfile but -1 expected!
  warn(
Traceback (most recent call last):
  File "/home/hieunguyen/miniconda3/envs/gctree/bin/gctree", line 7, in <module>
    sys.exit(main())
  File "/home/hieunguyen/miniconda3/envs/gctree/lib/python3.9/site-packages/gctree/cli.py", line 784, in main
    args.func(args)
  File "/home/hieunguyen/miniconda3/envs/gctree/lib/python3.9/site-packages/gctree/cli.py", line 169, in infer
    pp.parse_outfile(args.infiles[0], args.infiles[1], args.root)
  File "/home/hieunguyen/miniconda3/envs/gctree/lib/python3.9/site-packages/gctree/phylip_parse.py", line 193, in parse_outfile
    raise RuntimeError(f"No trees found in '{outfile}'")
RuntimeError: No trees found in './output/GCTree/IGHV1-47-01_IGHJ4-01_42_1/outfile'
mv: cannot stat 'IGHV1-47-01_IGHJ4-01_42_1*': No such file or directory
```
