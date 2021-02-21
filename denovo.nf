#!/usr/bin/env nextflow

if( params.help ) {
	return help()
}

reference_db = file(params.reference)
threads = params.threads
genome_size=params.genome_size
hybrid = file(params.hybrid)
meta_flag = params.metagenomic
long_only_flag = params.nanopore


Channel
    .fromFilePairs( params.reads, size: long_only_flag ? 1 : 2 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { reads }


process BwaIndexReference {
	input:
		path db_dir from reference_db
	output:
		file '*' into (host_reference_short, host_reference_long, host_reference_hybrid)

	"""
	bwa index $db_dir
	"""
}


process BwaHostDepletionShort {
    tag {samplename}

    input:
        set samplename, file(forward), file(reverse) from reads
        file(index) from host_reference_short
		file(reference_db)
	output:
		set samplename, file("${samplename}_short_R1.fastq"), file("${samplename}_short_R2.fastq") into (short_depleted)

	when:
		!long_only_flag

    """
	bwa mem -t $threads $reference_db $forward $reverse > ${samplename}.sam
	deplete_host_from_sam.py ${samplename}.sam -r1 ${samplename}_short_R1.fastq -r2 ${samplename}_short_R2.fastq
    """
}


process BwaHostDepletionLongExclusive {
    tag {samplename}

    input:
        set samplename, file(long) from reads
        file(index) from host_reference_long
		file(reference_db)
	output:
		set samplename, file("${samplename}_long.fastq") into (long_depleted)

	when:
		long_only_flag

    """
	bwa mem -t $threads $reference_db $long > ${samplename}.sam
	deplete_host_from_sam.py ${samplename}.sam -r1 ${samplename}_long.fastq
    """
}


process BwaHostDepletionLongHybrid {

    input:
        file(hybrid)
        file(index) from host_reference_hybrid
		file(reference_db)
	output:
		file("long.fastq") into (long_hybrid_depleted)

	when:
		hybrid && !long_only_flag

    """
	bwa mem -t $threads $reference_db $hybrid > long.sam
	deplete_host_from_sam.py long.sam -r1 long.fastq
    """
}


process SpadesAssemblyShort {
	tag {samplename}

	publishDir "${params.output}/SpadesAssemblies", mode: "copy"

	input:
		set samplename, file(forward), file(reverse) from short_depleted
	output:
		file("${samplename}_spades_assembly.fasta") into spades_out_short

	when:
		!long_only_flag

	script:
	if( meta_flag )
		"""
		spades.py -t $threads --tmp /tmp/spades_tmp -1 $forward -2 $reverse -o spades_output --meta
		mv spades_output/scaffolds.fasta ${samplename}_spades_assembly.fasta
		"""
	else
		"""
		spades -t $threads --tmp /tmp/spades_tmp -1 $forward -2 $reverse -o spades_output
		mv spades_output/scaffolds.fasta ${samplename}_spades_assembly.fasta
		"""
}


process SpadesAssemblyHybrid {
	tag {samplename}

	publishDir "${params.output}/SpadesAssemblies", mode: "copy"

	input:
		set samplename, file(forward), file(reverse) from short_depleted
		file(long) from long_hybrid_depleted
	output:
		file("${samplename}_spades_assembly.fasta") into spades_out_hybrid

	when:
		hybrid && !long_only_flag

	"""
	spades.py -t $threads --tmp /tmp/spades_tmp -1 $forward -2 $reverse --nanopore $long -o spades_output
	mv spades_output/scaffolds.fasta ${samplename}_spades_assembly.fasta
	"""
}


process FlyeAssemblyLong {
	tag {samplename}

	publishDir "${params.output}/FlyeAssemblies", mode: "copy"

	input:
		set samplename, file(long) from long_depleted
	output:
		file("${samplename}_flye_assembly.fasta") into flye_out
		file("${samplename}_flye_assembly_info.txt")

	when:
		long_only_flag

	script:
	if( meta_flag )
		"""
		flye --nano-raw $long -g ${genome_size} -t $threads -o flye_output --meta --resume
		"""
	else
		"""
		flye --nano-raw $long -g ${genome_size} -t $threads -o flye_output --resume
		mv flye_output/assembly.fasta ${samplename}_flye_assembly.fasta
		mv flye_output/assembly_info.txt ${samplename}_flye_assembly_info.txt
		"""
}


process GenomeEvaluationSpades {
	tag {samplename}

	publishDir "${params.output}/AssemblyStatistics", mode: "copy"

	input:
		file(spades_fasta) from spades_out_short.mix(spades_out_hybrid)
		file(reference_db)
	output:
		file "quast_output"

	when:
		!long_only_flag

	"""
	quast.py -t $threads -r $reference_db -o quast_output $spades_fasta
	"""
}


process GenomeEvaluationFlye {
	tag {samplename}

	publishDir "${params.output}/AssemblyStatistics", mode: "copy"

	input:
		file(flye_fasta) from flye_output
		file(reference_db)
	output:
		file "quast_output"

	when:
		long_only_flag

	"""
	quast.py -t $threads -r $reference_db -o quast_output $flye_fasta
	"""
}


def help() {
    println ""
    println "Program: denovo.nf"
    println "Developer: Steven Lakin, USDA APHIS"
    println "Description: De novo assembly pipeline for short paired-end, hybrid short/long, or Nanopore long data."
    println ""
    println "Usage: denovo.nf [options]"
    println "Example: denovo.nf"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to paired-end FASTQ input or concatenated Nanopore FASTQ input for a single sample"
    println "    --hybrid        STR      optional path to concatenated Nanopore FASTQ input for hybrid assembly"
    println "    --output        STR      path to output directory"
    println "    --reference     STR      path to pre-made BLAST database to query against"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of process threads, default 1 (max thread use = maxForks * threads)"
    println "    --metagenomic   FLAG     use metaSPAdes or Flye for metagenomic sample (cannot be used with hybrid assembly)"
    println "    --nanopore      FLAG     the input reads are concatenated Nanopore FASTQ data only (uses Flye assembler)"
    println "    --genome_size   STR      reference genome size (format [float][m/g] for mega/giga basepairs"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}
