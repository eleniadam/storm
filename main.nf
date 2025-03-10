#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.inputControl = null
params.inputTumor = null
params.inputDuplicons = null
params.inputPromoter = null
params.inputBand = null
params.inputSCopy = null

params.onOmics = null

params.samToolsContainer = null
params.samToolsCpu = null
params.samToolsMemory = null

params.bamToolsContainer = null
params.bamToolsCpu = null
params.bamToolsMemory = null

params.seqtkContainer = null
params.seqtkCpu = null
params.seqtkMemory = null

params.mrsfastContainer = null
params.mrsfastCpu = null
params.mrsfastMemory = null

params.telomhunterContainer = null
params.telomhunterCpu = null
params.telomhunterMemory = null

params.pyIdsContainer = null
params.pyIdsCpu = null
params.pyIdsMemory = null

params.blatContainer = null
params.blatCpu = null
params.blatMemory = null

publish_dir = { params.onOmics ? "/mnt/workflow/pubdir" : "." }
publish_dir_local = "."

process getIndexControl {
    container params.samToolsContainer
    cpus params.samToolsCpu
    memory params.samToolsMemory
    publishDir publish_dir

    input:
    path bamPath
    
    output:
    path "${bamPath}.bai"

    script:
    """
    samtools index -b ${bamPath} -@ 16
    """
}

process getIndexTumor {
    container params.samToolsContainer
    cpus params.samToolsCpu
    memory params.samToolsMemory
    publishDir publish_dir
    
    input:
    path bamPath
    
    output:
    path "${bamPath}.bai"

    script:
    """
    samtools index -b ${bamPath} -@ 16
    """
}

process convertToFastqControl {
    container params.bamToolsContainer
    cpus params.bamToolsCpu
    memory params.bamToolsMemory
    publishDir publish_dir
    
    input:
    path bamPath
    path indexPath

    output:
    path "${bamPath}.fastq"

    script:
    """
    bamtools convert -in ${bamPath} -format fastq -out ${bamPath}.fastq
    """
}

process convertToFastqTumor {
    container params.bamToolsContainer
    cpus params.bamToolsCpu
    memory params.bamToolsMemory
    publishDir publish_dir
    
    input:
    path bamPath
    path indexPath

    output:
    path "${bamPath}.fastq"

    script:
    """
    bamtools convert -in ${bamPath} -format fastq -out ${bamPath}.fastq
    """
}

process filterFastqControl {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir
    
    input:
    path fastqPath

    output:
    path "${fastqPath}_gte150bp.fastq"

    script:
    """
    seqtk seq -L 150 ${fastqPath} > ${fastqPath}_gte150bp.fastq
    """
}

process filterFastqTumor {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir
    
    input:
    path fastqPath

    output:
    path "${fastqPath}_gte150bp.fastq"

    script:
    """
    seqtk seq -L 150 ${fastqPath} > ${fastqPath}_gte150bp.fastq
    """
}

process getIndexDuplicons {
    container params.mrsfastContainer
    cpus params.mrsfastCpu
    memory params.mrsfastMemory
    publishDir publish_dir

    input:
    path fastaPath

    output:
    path "${fastaPath}.index"

    script:
    """
    mrsfast --index ${fastaPath}
    """
}

process getDupliconMappingsControl {
    container params.mrsfastContainer
    cpus params.mrsfastCpu
    memory params.mrsfastMemory
    publishDir publish_dir

    input:
    path fastaPath
    path fastqPath
    path indexPath

    output:
    path "${fastqPath}_smappings.sam"

    script:
    """
    mrsfast --search ${fastaPath} --seq ${fastqPath} -e 10 --threads 4 -o ${fastqPath}_smappings.sam --disable-nohits
    """
}

process getDupliconMappingsTumor {
    container params.mrsfastContainer
    cpus params.mrsfastCpu
    memory params.mrsfastMemory
    publishDir publish_dir

    input:
    path fastaPath
    path fastqPath
    path indexPath

    output:
    path "${fastqPath}_smappings.sam"

    script:
    """
    mrsfast --search ${fastaPath} --seq ${fastqPath} -e 10 --threads 4 -o ${fastqPath}_smappings.sam --disable-nohits
    """
}

process getTelomericMappings {
    container params.telomhunterContainer
    cpus params.telomhunterCpu
    memory params.telomhunterMemory
    publishDir publish_dir

    input:
    path bamPathControl
    path bamPathTumor
    path bandPath
    path indexControl
    path indexTumor

    output:
    path "out_telomerehunter"

    script:
    """
    telomerehunter -ibc ${bamPathControl} -ibt ${bamPathTumor} -b ${bandPath} -o out_telomerehunter -p teloms -pl
    """
}

process getDupliconReadsControl {
    container params.samToolsContainer
    cpus params.samToolsCpu
    memory params.samToolsMemory
    publishDir publish_dir
    
    input:
    path samPath
    
    output:
    path "reads_subtelom_control.fastq"

    script:
    """
    samtools fastq ${samPath} > reads_subtelom_control.fastq
    """
}

process getDupliconReadsTumor {
    container params.samToolsContainer
    cpus params.samToolsCpu
    memory params.samToolsMemory
    publishDir publish_dir
    
    input:
    path samPath
    
    output:
    path "reads_subtelom_tumor.fastq"

    script:
    """
    samtools fastq ${samPath} > reads_subtelom_tumor.fastq
    """
}

process getDupliconMateIDsControl {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path fastqPath

    output:
    path "${fastqPath}_smateIDs.lst"

    script:
    """
    python3 /extract_read_ids.py mates fastq ${fastqPath} ${fastqPath}_smateIDs.lst
    """
}

process getDupliconMateIDsTumor {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path fastqPath

    output:
    path "${fastqPath}_smateIDs.lst"

    script:
    """
    python3 /extract_read_ids.py mates fastq ${fastqPath} ${fastqPath}_smateIDs.lst
    """
}

process getDupliconMatesControl {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "mates_subtelom_control.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > mates_subtelom_control.fastq
    """
}

process getDupliconMatesTumor {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "mates_subtelom_tumor.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > mates_subtelom_tumor.fastq
    """
}

process getTelomReadsControl {
    container params.bamToolsContainer
    cpus params.bamToolsCpu
    memory params.bamToolsMemory
    publishDir publish_dir
    
    input:
    path folderPath

    output:
    path "reads_telom_control.fastq"

    script:
    """
    bamtools convert -in ${folderPath}/teloms/control_TelomerCnt_teloms/teloms_filtered.bam -format fastq -out reads_telom_control.fastq
    """
}

process getTelomReadsTumor {
    container params.bamToolsContainer
    cpus params.bamToolsCpu
    memory params.bamToolsMemory
    publishDir publish_dir
    
    input:
    path folderPath

    output:
    path "reads_telom_tumor.fastq"

    script:
    """
    bamtools convert -in ${folderPath}/teloms/tumor_TelomerCnt_teloms/teloms_filtered.bam -format fastq -out reads_telom_tumor.fastq
    """
}

process getTelomMateIDsControl {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path fastqPath

    output:
    path "${fastqPath}_tmateIDs.lst"

    script:
    """
    python3 /extract_read_ids.py mates fastq ${fastqPath} ${fastqPath}_tmateIDs.lst
    """
}

process getTelomMateIDsTumor {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path fastqPath

    output:
    path "${fastqPath}_tmateIDs.lst"

    script:
    """
    python3 /extract_read_ids.py mates fastq ${fastqPath} ${fastqPath}_tmateIDs.lst
    """
}

process getTelomMatesControl {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "mates_telom_control.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > mates_telom_control.fastq
    """
}

process getTelomMatesTumor {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "mates_telom_tumor.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > mates_telom_tumor.fastq
    """
}

process convertToFastaControl {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath

    output:
    path "${fastqPath}.fasta"

    script:
    """
    seqtk seq -a ${fastqPath} > ${fastqPath}.fasta
    """
}

process convertToFastaTumor {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath

    output:
    path "${fastqPath}.fasta"

    script:
    """
    seqtk seq -a ${fastqPath} > ${fastqPath}.fasta
    """
}

process getPromoterMappingsControl {
    container params.blatContainer
    cpus params.blatCpu
    memory params.blatMemory
    publishDir publish_dir

    input:
    path fastaPromoterPath
    path fastaPath

    output:
    path "${fastaPath}_prommappings.out"

    script:
    """
    blat ${fastaPromoterPath} ${fastaPath} ${fastaPath}_prommappings.out -t=dna -q=dna -out=blast8
    """
}

process getPromoterMappingsTumor {
    container params.blatContainer
    cpus params.blatCpu
    memory params.blatMemory
    publishDir publish_dir

    input:
    path fastaPromoterPath
    path fastaPath

    output:
    path "${fastaPath}_prommappings.out"

    script:
    """
    blat ${fastaPromoterPath} ${fastaPath} ${fastaPath}_prommappings.out -t=dna -q=dna -out=blast8
    """
}

process getPromoterReadIDsControl {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path txtPath

    output:
    path "${txtPath}_preadIDs.lst"

    script:
    """
    python3 /extract_read_ids.py reads blat ${txtPath} ${txtPath}_preadIDs.lst
    """
}

process getPromoterReadIDsTumor {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path txtPath

    output:
    path "${txtPath}_preadIDs.lst"

    script:
    """
    python3 /extract_read_ids.py reads blat ${txtPath} ${txtPath}_preadIDs.lst
    """
}

process getPromoterMateIDsControl {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path txtPath

    output:
    path "${txtPath}_pmateIDs.lst"

    script:
    """
    python3 /extract_read_ids.py mates blat ${txtPath} ${txtPath}_pmateIDs.lst
    """
}

process getPromoterMateIDsTumor {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path txtPath

    output:
    path "${txtPath}_pmateIDs.lst"

    script:
    """
    python3 /extract_read_ids.py mates blat ${txtPath} ${txtPath}_pmateIDs.lst
    """
}

process getPromoterReadsControl {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "reads_prom_control.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > reads_prom_control.fastq
    """
}

process getPromoterReadsTumor {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "reads_prom_tumor.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > reads_prom_tumor.fastq
    """
}

process getPromoterMatesControl {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "mates_prom_control.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > mates_prom_control.fastq
    """
}

process getPromoterMatesTumor {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "mates_prom_tumor.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > mates_prom_tumor.fastq
    """
}

process getIndexSCopy {
    container params.mrsfastContainer
    cpus params.mrsfastCpu
    memory params.mrsfastMemory
    publishDir publish_dir

    input:
    path fastaPath

    output:
    path "${fastaPath}.index"

    script:
    """
    mrsfast --index ${fastaPath}
    """
}

process getSCopyMappingsControl {
    container params.mrsfastContainer
    cpus params.mrsfastCpu
    memory params.mrsfastMemory
    publishDir publish_dir

    input:
    path fastaPath
    path fastqPath
    path indexPath

    output:
    path "${fastqPath}_scpymappings.sam"

    script:
    """
    mrsfast --search ${fastaPath} --seq ${fastqPath} -e 10 --threads 4 -o ${fastqPath}_scpymappings.sam --disable-nohits
    """
}

process getSCopyMappingsTumor {
    container params.mrsfastContainer
    cpus params.mrsfastCpu
    memory params.mrsfastMemory
    publishDir publish_dir

    input:
    path fastaPath
    path fastqPath
    path indexPath

    output:
    path "${fastqPath}_scpymappings.sam"

    script:
    """
    mrsfast --search ${fastaPath} --seq ${fastqPath} -e 10 --threads 4 -o ${fastqPath}_scpymappings.sam --disable-nohits
    """
}

process getSCopyReadsControl {
    container params.samToolsContainer
    cpus params.samToolsCpu
    memory params.samToolsMemory
    publishDir publish_dir
    
    input:
    path samPath
    
    output:
    path "reads_scopy_control.fastq"

    script:
    """
    samtools fastq ${samPath} > reads_scopy_control.fastq
    """
}

process getSCopyReadsTumor {
    container params.samToolsContainer
    cpus params.samToolsCpu
    memory params.samToolsMemory
    publishDir publish_dir
    
    input:
    path samPath
    
    output:
    path "reads_scopy_tumor.fastq"

    script:
    """
    samtools fastq ${samPath} > reads_scopy_tumor.fastq
    """
}

process getSCopyMateIDsControl {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path fastqPath

    output:
    path "${fastqPath}_scpymateIDs.lst"

    script:
    """
    python3 /extract_read_ids.py mates fastq ${fastqPath} ${fastqPath}_scpymateIDs.lst
    """
}

process getSCopyMateIDsTumor {
    container params.pyIdsContainer
    cpus params.pyIdsCpu
    memory params.pyIdsMemory
    publishDir publish_dir_local

    input:
    path fastqPath

    output:
    path "${fastqPath}_scpymateIDs.lst"

    script:
    """
    python3 /extract_read_ids.py mates fastq ${fastqPath} ${fastqPath}_scpymateIDs.lst
    """
}

process getSCopyMatesControl {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "mates_scopy_control.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > mates_scopy_control.fastq
    """
}

process getSCopyMatesTumor {
    container params.seqtkContainer
    cpus params.seqtkCpu
    memory params.seqtkMemory
    publishDir publish_dir

    input:
    path fastqPath
    path listPath

    output:
    path "mates_scopy_tumor.fastq"

    script:
    """
    seqtk subseq ${fastqPath} ${listPath} > mates_scopy_tumor.fastq
    """
}

workflow {
    outIndexDuplicons = getIndexDuplicons(params.inputDuplicons)
    
    outIndexControl = getIndexControl(params.inputControl)
    outFastqControl = convertToFastqControl(params.inputControl, outIndexControl)
    outFilteredFastqControl = filterFastqControl(outFastqControl)
    outDupliconMappingsControl = getDupliconMappingsControl(params.inputDuplicons, outFilteredFastqControl, outIndexDuplicons)
    outDupliconReadsControl = getDupliconReadsControl(outDupliconMappingsControl)
    outDupliconMateIDsControl = getDupliconMateIDsControl(outDupliconReadsControl)
    getDupliconMatesControl(outFilteredFastqControl, outDupliconMateIDsControl)
    
    outIndexTumor = getIndexTumor(params.inputTumor)
    outFastqTumor = convertToFastqTumor(params.inputTumor, outIndexTumor)
    outFilteredFastqTumor = filterFastqTumor(outFastqTumor)
    outDupliconMappingsTumor = getDupliconMappingsTumor(params.inputDuplicons, outFilteredFastqTumor, outIndexDuplicons)
    outDupliconReadsTumor = getDupliconReadsTumor(outDupliconMappingsTumor)
    outDupliconMateIDsTumor = getDupliconMateIDsTumor(outDupliconReadsTumor)
    getDupliconMatesTumor(outFilteredFastqTumor, outDupliconMateIDsTumor)
    
    outMappingsTelom = getTelomericMappings(params.inputControl, params.inputTumor, params.inputBand, outIndexControl, outIndexTumor)
    outTelomReadsControl = getTelomReadsControl(outMappingsTelom)
    outTelomMateIDsControl = getTelomMateIDsControl(outTelomReadsControl)
    getTelomMatesControl(outFilteredFastqControl, outTelomMateIDsControl)
    outTelomReadsTumor = getTelomReadsTumor(outMappingsTelom)
    outTelomMateIDsTumor = getTelomMateIDsTumor(outTelomReadsTumor)
    getTelomMatesTumor(outFilteredFastqTumor, outTelomMateIDsTumor)

    outFilteredFastaControl = convertToFastaControl(outFilteredFastqControl)
    outMappingsPromControl = getPromoterMappingsControl(params.inputPromoter, outFilteredFastaControl)
    outPromReadIDsControl = getPromoterReadIDsControl(outMappingsPromControl)
    getPromoterReadsControl(outFilteredFastqControl, outPromReadIDsControl)
    outPromMateIDsControl = getPromoterMateIDsControl(outMappingsPromControl)
    getPromoterMatesControl(outFilteredFastqControl, outPromMateIDsControl)
    outFilteredFastaTumor = convertToFastaTumor(outFilteredFastqTumor)
    outMappingsPromTumor = getPromoterMappingsTumor(params.inputPromoter, outFilteredFastaTumor)
    outPromReadIDsTumor = getPromoterReadIDsTumor(outMappingsPromTumor)
    getPromoterReadsTumor(outFilteredFastqTumor, outPromReadIDsTumor)
    outPromMateIDsTumor = getPromoterMateIDsTumor(outMappingsPromTumor)
    getPromoterMatesTumor(outFilteredFastqTumor, outPromMateIDsTumor)

    outIndexSCopy = getIndexSCopy(params.inputSCopy)
    outScopyMappingsControl = getSCopyMappingsControl(params.inputSCopy, outFilteredFastqControl, outIndexSCopy)
    outSCopyReadsControl = getSCopyReadsControl(outScopyMappingsControl)
    outSCopyMateIDsControl = getSCopyMateIDsControl(outSCopyReadsControl)
    getSCopyMatesControl(outFilteredFastqControl, outSCopyMateIDsControl)
    outScopyMappingsTumor = getSCopyMappingsTumor(params.inputSCopy, outFilteredFastqTumor, outIndexSCopy)
    outSCopyReadsTumor = getSCopyReadsTumor(outScopyMappingsTumor)
    outSCopyMateIDsTumor = getSCopyMateIDsTumor(outSCopyReadsTumor)
    getSCopyMatesTumor(outFilteredFastqTumor, outSCopyMateIDsTumor)
}


