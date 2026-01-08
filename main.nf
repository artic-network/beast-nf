#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * BEAST X Phylogenetic Analysis Pipeline
 * 
 * This pipeline performs:
 * 1. BEAST XML generation from aligned FASTA
 * 2. BEAST X phylogenetic analysis
 * 3. TreeAnnotator MCC tree summary
 * 4. Time tree visualization with Baltic
 */

// Parameters
params.input = null
params.template = null
params.outdir = "results"
params.prefix = "beast_analysis"
params.chain_length = 10000000
params.log_every = 1000
params.screen_every = 10000
params.burnin = 10

// Validate input
if (!params.input) {
    error "Please provide an input FASTA file with --input"
}
if (!params.template) {
    error "Please provide a BEAST XML template file with --template"
}

log.info """
    BEAST X ANALYSIS PIPELINE
    =========================
    input        : ${params.input}
    template     : ${params.template}
    outdir       : ${params.outdir}
    prefix       : ${params.prefix}
    chain_length : ${params.chain_length}
    burnin       : ${params.burnin}%
    """
    .stripIndent()

/*
 * PROCESS: Generate BEAST XML
 */
process GENERATE_XML {
    publishDir "${params.outdir}/xml", mode: 'copy'
    
    input:
    path fasta
    path template
    
    output:
    path "${params.prefix}.xml", emit: xml
    
    script:
    """
    beastgen \\
        -D chain_length=${params.chain_length} \\
        -D log_every=${params.log_every} \\
        -D screen_every=${params.screen_every} \\
        ${template} \\
        ${fasta} \\
        ${params.prefix}.xml
    """
}

/*
 * PROCESS: Run BEAST X
 */
process RUN_BEAST {
    publishDir "${params.outdir}/beast", mode: 'copy'
    
    input:
    path xml
    
    output:
    path "${params.prefix}.log", emit: log
    path "${params.prefix}.trees", emit: trees
    path "${params.prefix}.*.log", emit: all_logs
    
    script:
    """
    beast ${xml}
    """
}

/*
 * PROCESS: TreeAnnotator
 */
process TREE_ANNOTATOR {
    publishDir "${params.outdir}/trees", mode: 'copy'
    
    input:
    path trees
    
    output:
    path "${params.prefix}.mcc.tree", emit: mcc_tree
    
    script:
    """
    treeannotator \\
        --burnin ${params.burnin} \\
        --type hipstr \\
        ${trees} \\
        ${params.prefix}.mcc.tree
    """
}

/*
 * PROCESS: Visualize with Baltic
 */
process VISUALIZE_TREE {
    publishDir "${params.outdir}/figures", mode: 'copy'
    
    input:
    path mcc_tree
    
    output:
    path "${params.prefix}_timetree.png", emit: png
    path "${params.prefix}_timetree.svg", emit: svg
    
    script:
    """
    visualize_tree.py \\
        --input ${mcc_tree} \\
        --output ${params.prefix}_timetree \\
        --format png,svg
    """
}

/*
 * PROCESS: Generate HTML Report
 */
process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path fasta
    path template
    path log
    path svg
    val runtime
    
    output:
    path "${params.prefix}_report.html", emit: report
    
    script:
    """
    # Run loganalyser
    loganalyser -burnin ${params.burnin} ${log} > loganalyser_output.txt
    
    # Generate HTML report
    generate_report.py \\
        --fasta ${fasta} \\
        --template ${template} \\
        --log ${log} \\
        --loganalyser loganalyser_output.txt \\
        --svg ${svg} \\
        --output ${params.prefix}_report.html \\
        --chain-length ${params.chain_length} \\
        --log-every ${params.log_every} \\
        --burnin ${params.burnin} \\
        --runtime '${runtime}'
    """
}

/*
 * WORKFLOW
 */
workflow {
    // Create channel from input FASTA file
    fasta_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Create channel from template file
    template_ch = Channel.fromPath(params.template, checkIfExists: true)
    
    // Generate BEAST XML
    GENERATE_XML(fasta_ch, template_ch)
    
    // Run BEAST
    RUN_BEAST(GENERATE_XML.out.xml)
    
    // Run TreeAnnotator
    TREE_ANNOTATOR(RUN_BEAST.out.trees)
    
    // Visualize tree
    VISUALIZE_TREE(TREE_ANNOTATOR.out.mcc_tree)
    
    // Generate report
    GENERATE_REPORT(
        fasta_ch,
        template_ch,
        RUN_BEAST.out.log,
        VISUALIZE_TREE.out.svg,
        workflow.duration
    )
}

workflow.onComplete {
    log.info """
        Pipeline completed!
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration  : ${workflow.duration}
        Results   : ${params.outdir}
        """
        .stripIndent()
}
