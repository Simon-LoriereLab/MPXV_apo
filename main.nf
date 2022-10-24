nextflow.enable.dsl=2

params.reference="data/reference.fasta"
params.results="results"
params.data='data/batch1/*.fas'

process cutReference {
    label 'goalign'
    publishDir "${params.results}/", mode: 'copy'

    input:
    path ref

    output:
    path "${ref.baseName}_cut.fasta"

    script:
    """
    goalign subseq -s 168239 -l 248 -i $ref -o ${ref.baseName}_cut.fasta
    """
}

process rename {
    label 'goalign'
    tag "$fa"
    
    input:
    path fa
    
    output:
    path "*_rename.fasta"
    
    script:
    """
    goalign rename -e ' ' -b '_' -i $fa -o ${fa.baseName.replaceAll("\\s","_")}_rename.fasta
    """
}

process minimap {
    label 'map'
    tag "$fa"

    input:
    path reference
    path fa

    output:
    path "*.sam"

    script:
    """
    minimap2 -a -x splice --sam-hit-only --secondary=no --score-N=0 ${reference} ${fa} -o ${fa.baseName}.sam
    """
}

process toFasta {
    label 'gofasta'
    tag "$sam"

    input:
    path sam

    output:
    path "*.fasta"

    script:
    """
    gofasta sam toMultiAlign -s ${sam} --start 168240 --end 168487 -o ${sam.baseName}_al.fasta
    # --pad
    """
}

process cleanseq {
    label 'goalign'
    tag "${fa}"
    publishDir "${params.results}", mode: 'copy'

    input:
    path fa

    output:
    path "${fa.baseName}_clean.fasta"

    script:
    """
    goalign clean seqs --char N -c 0.7 -i ${fa} -o ${fa.baseName}_clean.fasta
    """
}

process countMutations {
    label 'python'
    tag "${align}"

    publishDir "${params.results}", mode: 'copy'
    
    input:
    path ref
    path align

    output:
    path "mutations.txt"

    script:
    """
    countMutations.py -i ${align} -r ${ref} > mutations.txt
    """
}

workflow {
    ref = file(params.reference)
    fa = Channel.fromPath(params.data)

    refcut=cutReference(ref)
    ren=rename(fa)
    sam=minimap(ref,ren)
    alfa=toFasta(sam)
    alfaall=alfa.collectFile(name:"align.fasta", storeDir: "${params.results}/")
    cleanalfaall=cleanseq(alfaall)
    countMutations(refcut,cleanalfaall)
}