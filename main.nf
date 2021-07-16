nextflow.enable.dsl = 2

helpMessage = """
Use short reads to polish an assembly.

USAGE: nextflow run WarrenLab/shortread-polish-nf
    --assembly [fasta]
    --sra [accession] OR --fastq [path]
"""

params.maxCoverage = 500

process faidx {
    input:
    path assembly

    output:
    path "${assembly}.fai", emit: fai

    """
    samtools faidx $assembly
    """
}

process align {
    publishDir 'bams'

    input:
    path assembly
    tuple(val(sra_id), file(reads))

    output:
    file "${sra_id}.bam"

    """
    minimap2 -ax sr -t ${task.cpus} ${assembly} ${reads} | \
        samtools view -bh - | samtools fixmate -m - - | \
        samtools sort -T . -@ ${task.cpus} -m 2G - | \
        samtools markdup - ${sra_id}.bam
    """
}

process mergeBams {
    input:
    file bams

    output:
    file("merged.bam")

    """
    bams="${bams}"
    words=( \$bams )
    if [[ "1" == "\${#words[@]}" ]]; then
        cp ${bams} merged.bam
    else
        samtools merge merged.bam ${bams}
    fi
    """
}


process freebayes {
    input:
    path assembly
    path "merged.bam"
    each region

    output:
    file "${region.name}.bcf"

    """
    samtools index merged.bam
    freebayes --bam merged.bam --region ${region.name}:1-${region.length} \
        --skip-coverage 500 -f ${assembly} \
        | bcftools view --no-version -Ob -o ${region.name}.bcf
    """
}

process concatAndConsensus {
    publishDir 'consensus'

    input:
    path assembly
    path "*.bcf"

    output:
    file "polished.fa"
    file "joined.bcf"
    file "report.txt"

    """
    bcftools concat -n *.bcf | bcftools view -Ou -e'type="ref"' \
        | bcftools norm -Ob -f $assembly -o joined.bcf
    bcftools index joined.bcf
    bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Hla -f $assembly joined.bcf > polished.fa

    snp=\$(bcftools view -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Ha joined.bcf | grep -c 'TYPE=snp')
    ins=\$(bcftools view -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Ha joined.bcf | grep -c 'TYPE=ins')
    del=\$(bcftools view -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Ha joined.bcf | grep -c 'TYPE=del')
    echo "\${snp} SNPs, \${ins} insertions, and \${del} deletions corrected." \
        > report.txt
    """
}

workflow {
    if (params.help) {
        print helpMessage
        exit 0
    }

    if (params.sra) {
        shortReads = Channel.fromSRA(params.sra)
    } else if (params.fastq) {
        shortReads = Channel.fromFilePairs(params.fastq)
    } else {
        println "ERROR: Need either --sra or --fastq option!"
        print helpMessage
    }

    assembly = Channel.fromPath(params.assembly)

    faidx(assembly)

    regions = faidx.out.fai
        .splitCsv(header: ['name', 'length', 'a', 'b', 'c'], sep: "\t")

    align(assembly, shortReads)
    mergeBams(align.out.collect())
    freebayes(assembly, mergeBams.out, regions)
    concatAndConsensus(assembly, freebayes.out.collect())
}
