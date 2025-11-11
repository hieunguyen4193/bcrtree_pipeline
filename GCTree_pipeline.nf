nextflow.enable.dsl=1

// path to fetch fasta file as inputs
params.samplesheet=""
// Validate sample sheet parameter
if (params.samplesheet == '') {
    log.error "Sample sheet not specified. Use --samplesheet to specify the path to your sample sheet."
    System.exit(1)
}

// path to save output all output files
params.output=""

// path to the python source use in deduplicating the sequences.
deduplicate_src=file(params.deduplicate_src)
// modify_tree_colors=file(params.modify_tree_colors)
// color_path=file(params.color_path)

// Create a channel that emits files from the specified input directory
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: ',', strip: true)
    .map { row -> 
        def sample_id = row.filename
        def path = row.path
        return [sample_id, file(path)]
    }
    .set { fastaFilesChannel }
    
process deduplicate { 
    cache "deep"; tag "$sample_id"
    publishDir "$params.output/$sample_id/01_deduplicate", mode: 'copy'
    errorStrategy 'terminate'
    maxRetries 1
    maxForks 20

    input:
        tuple sample_id, path(fasta) from fastaFilesChannel
        file(deduplicate_src)
    output: 
        tuple sample_id, file("${sample_id}*") into mkconfig_ch
        tuple sample_id, file(fasta) into orig_fasta_ch
        tuple sample_id, "${sample_id}.abundance.csv" into abundance_ch
        tuple sample_id, "${sample_id}.id_map.csv" into idmap_ch
        tuple sample_id, "${sample_id}.id_map_seq.csv" into idmap_seq_ch
        tuple sample_id, "${sample_id}.phylip" into phylip_ch

    shell:
    '''
    cat !{fasta} | sed "s/>.*|Abundance:\\([0-9]\\+\\)/>\\1/" > !{fasta}.modified.fa
    python !{deduplicate_src} \
        --input !{fasta}.modified.fa \
        --root GL \
        --frame 0 \
        --id_abundances \
        --output_name !{sample_id} \
        --output . 
    '''
}

process dnapars_and_inferring_gc_trees {
    cache "deep"; tag "$sample_id"
    publishDir "$params.output/$sample_id/02_dnapars", mode: 'copy'
    errorStrategy 'ignore'
    // errorStrategy 'terminate'
    maxRetries 1
    maxForks 20
    
    input:
        tuple sample_id, file("${sample_id}*") from mkconfig_ch
        tuple sample_id, "${sample_id}.abundance.csv" from abundance_ch
        tuple sample_id, "${sample_id}.phylip" from phylip_ch
    output:
        tuple sample_id, file("*") into modify_gctree_colors_ch

    shell:
    '''
    mkconfig --quick !{sample_id}.phylip dnapars > !{sample_id}_dnapars.cfg
    dnapars < !{sample_id}_dnapars.cfg > !{sample_id}_dnapars.log

    export QT_QPA_PLATFORM=offscreen
    export XDG_RUNTIME_DIR=/tmp/runtime-runner
    export MPLBACKEND=agg

    gctree infer --verbose --root GL --frame 1 --idlabel outfile !{sample_id}.abundance.csv
    '''
}