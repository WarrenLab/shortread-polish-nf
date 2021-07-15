nextflow.enable.dsl = 2

// necessary parameters:
// --assembly [fasta]: assembly to polish
// --sra [accession] OR --fastq [path]: short reads
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
    tuple(file("merged.bam"), val(region))

    output:
    file "${region.name}.bcf"

    """
    samtools index merged.bam
    freebayes --bam merged.bam --region ${region.name}:1-${region.length} \
        --skip-coverage 500 -f ${params.assembly} \
        | bcftools view --no-version -Ob -o ${region.name}.bcf
    """
}

process concat_and_consensus {
    publishDir 'consensus'

    input:
    file "*.bcf"

    output:
    file "polished.fa"
    file "joined.bcf"
    file "report.txt"

    """
    bcftools concat -n *.bcf | bcftools view -Ou -e'type="ref"' \
        | bcftools norm -Ob -f ${params.assembly} -o joined.bcf
    bcftools index joined.bcf
    bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' \
        -Hla -f ${params.assembly} joined.bcf > polished.fa

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
    if (params.sra) {
        shortReads = Channel.fromSRA(params.sra)
    } else if (params.fastq) {
        shortReads = Channel.fromFilePairs(params.fastq)
    } else {
        println("Need either --sra or --fastq option!")
    }
    
    assembly = Channel.fromPath(params.assembly)

    faidx(assembly)

    regions = faidx.out.fai
        .splitCsv(header: ['name', 'length', 'a', 'b', 'c'], sep: "\t")

    align(assembly, shortReads)
    mergeBams(align.out.collect())
    
    concatAndConsensus(freebayes(mergeBams.out.combine(regions)).out.collect())
}

