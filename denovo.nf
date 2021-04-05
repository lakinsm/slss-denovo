#!/usr/bin/env nextflow

if( params.help ) {
	return help()
}

diamond_db = file(params.diamond_db)
reference_db = file(params.reference)
host_reference = file(params.host)
threads = params.threads
hybrid = file(params.hybrid)
reverse = file(params.reverse)
meta_flag = params.metagenomic
long_only_flag = params.nanopore
my_scratch = params.scratch_dir


Channel
    .fromPath( params.forward )
    .map{ in -> tuple( in.getName(), file(in)) }
    .into{ forward_reads; nanopore_reads }


Channel
	.fromPath( params.reverse )
	.set{ reverse_reads }


Channel
	.fromPath( params.hybrid )
	.set{ hybrid_reads }


if( params.use_index ) {
	host_reference_short = Channel.fromPath("${params.host}.*")
	host_reference_long = Channel.fromPath("${params.host}.*")
	host_reference_hybrid = Channel.fromPath("${params.host}.*")
} else {
	process BwaIndexReference {
		input:
			path db_dir from host_reference
		output:
			file "${params.host}.*" into (host_reference_short, host_reference_long, host_reference_hybrid)

		script:
		if( db_dir.name != 'NONE_REF' )
			"""
			bwa index $db_dir
			"""
		else
			"""
			touch empty.fasta
			"""
	}
}



process BwaHostDepletionShort {
    tag {samplename}

    publishDir "${params.output}/HostDepletion", mode: "copy"

    input:
        set samplename, path(forward) from forward_reads
        file(reverse)
        file index from host_reference_short.collect()
		file opt_ref from host_reference
	output:
		set samplename, file("${samplename}_short_R1.fastq"), file("${samplename}_short_R2.fastq") into (short_depleted, short_hybrid_depleted)
		file("*.log")

	when:
		!long_only_flag

	script:
	if( opt_ref.name != 'NONE_REF' )
	    """
		bwa mem -t $threads $opt_ref $forward $reverse > ${samplename}.sam
		deplete_host_from_sam.py ${samplename}.sam -r1 ${samplename}_short_R1.fastq -r2 ${samplename}_short_R2.fastq
	    """
	else
		"""
		cp forward ${samplename}_short_R1.fastq
		cp reverse ${samplename}_short_R2.fastq
		echo "No depletion performed" > ${samplename}_depletion_stats.log
		"""
}


process BwaHostDepletionLongExclusive {
    tag {samplename}

    publishDir "${params.output}/HostDepletion", mode: "copy"

    input:
        set samplename, path(reads) from nanopore_reads
        file index from host_reference_long.collect()
		file opt_ref from host_reference
	output:
		set samplename, file("${samplename}_long.fastq") into (long_depleted)
		file("*.log")

	when:
		long_only_flag

	script:
	if( opt_ref.name != 'NONE_REF' )
	    """
		bwa mem -t $threads $opt_ref $reads > ${samplename}.sam
		deplete_host_from_sam.py ${samplename}.sam -r1 ${samplename}_long.fastq
	    """
	else
		"""
		cp $reads ${samplename}_long.fastq
		echo "No depletion performed" > ${samplename}_depletion_stats.log
		"""
}


process BwaHostDepletionLongHybrid {
	tag {samplename}

	publishDir "${params.output}/HostDepletion", mode: "copy"

    input:
        file hybrid_str from hybrid
        file index from host_reference_hybrid.collect()
		file opt_ref from host_reference
	output:
		file("long.fastq") into (long_hybrid_depleted)
		file("*.log")

	when:
		(hybrid_str.name != 'NONE_HYBRID') && !long_only_flag

	script:
	if( opt_ref.name != 'NONE_REF' )
	    """
		bwa mem -t $threads $opt_ref $hybrid > long.sam
		deplete_host_from_sam.py long.sam -r1 long.fastq
	    """
	else
		"""
		cp $hybrid long.fastq
		echo "No depletion performed" > long_depletion_stats.log
		"""
}


process SpadesAssemblyShort {
	tag {samplename}

	publishDir "${params.output}/SpadesAssemblies", mode: "copy"

	input:
		set samplename, file(forward), file(reverse) from short_depleted
		file hybrid_str from hybrid
	output:
		file("${samplename}_spades_assembly.fasta") into (quast_spades_out_short, dmnd_spades_out_short)

	when:
		(hybrid_str.name == 'NONE_HYBRID') && !long_only_flag

	script:
	if( meta_flag )
		"""
		spades.py -t $threads --tmp $my_scratch -1 $forward -2 $reverse -o spades_output --meta
		mv spades_output/scaffolds.fasta ${samplename}_spades_assembly.fasta
		"""
	else
		"""
		spades.py -t $threads --tmp $my_scratch -1 $forward -2 $reverse -o spades_output
		mv spades_output/scaffolds.fasta ${samplename}_spades_assembly.fasta
		"""
}


process SpadesAssemblyHybrid {
	tag {samplename}

	publishDir "${params.output}/SpadesAssemblies", mode: "copy"

	input:
		set samplename, file(forward), file(reverse) from short_hybrid_depleted
		file(long_reads) from long_hybrid_depleted
		file hybrid_str from hybrid
	output:
		file("${samplename}_spades_assembly.fasta") into (quast_spades_out_hybrid, dmnd_spades_out_hybrid)

	when:
		(hybrid_str.name != 'NONE_HYBRID') && !long_only_flag

	script:
	if( meta_flag )
		"""
		spades.py -t $threads --tmp $my_scratch -1 $forward -2 $reverse --nanopore $long_reads -o spades_output --meta
		mv spades_output/scaffolds.fasta ${samplename}_spades_assembly.fasta
		"""
	else
		"""
		spades.py -t $threads --tmp $my_scratch -1 $forward -2 $reverse --nanopore $long -o spades_output
		mv spades_output/scaffolds.fasta ${samplename}_spades_assembly.fasta
		"""
}


process FlyeAssemblyLong {
	tag {samplename}

	publishDir "${params.output}/FlyeAssemblies", mode: "copy"

	input:
		set samplename, file(reads) from long_depleted
	output:
		file("${samplename}_flye_assembly.fasta") into (quast_flye_out, dmnd_flye_out)
		file("${samplename}_flye_assembly_info.txt")

	when:
		long_only_flag

	script:
	if( meta_flag )
		"""
		flye --nano-raw $reads -t $threads -o flye_output --meta
		mv flye_output/assembly.fasta ${samplename}_flye_assembly.fasta
		mv flye_output/assembly_info.txt ${samplename}_flye_assembly_info.txt
		"""
	else
		"""
		flye --nano-raw $reads -t $threads -o flye_output
		mv flye_output/assembly.fasta ${samplename}_flye_assembly.fasta
		mv flye_output/assembly_info.txt ${samplename}_flye_assembly_info.txt
		"""
}


process GenomeEvaluationSpades {
	publishDir "${params.output}/AssemblyStatistics", mode: "copy"

	input:
		file(spades_fasta) from quast_spades_out_short.mix(quast_spades_out_hybrid)
		file opt_ref from reference_db
	output:
		file "quast_output"

	when:
		!long_only_flag

	script:
	if( opt_ref.name != 'NONE_REF' )
		"""
		quast.py -t $threads -r $opt_ref -o quast_output $spades_fasta
		"""
	else
		"""
		quast.py -t $threads -o quast_output $spades_fasta
		"""
}


process GenomeEvaluationFlye {
	publishDir "${params.output}/AssemblyStatistics", mode: "copy"

	input:
		file(flye_fasta) from quast_flye_out
		file opt_ref from reference_db
	output:
		file "quast_output"

	when:
		long_only_flag

	script:
	if( opt_ref.name != 'NONE_REF' )
		"""
		quast.py -t $threads -r $opt_ref -o quast_output $flye_fasta
		"""
	else
		"""
		quast.py -t $threads -o quast_output $spades_fasta
		"""
}


process DiamondBlast {
	publishDir "${params.output}/DiamondBlast", mode: "copy"

	input:
		file(assembly_fasta) from dmnd_flye_out.mix(dmnd_spades_out_hybrid, dmnd_spades_out_short)
		file(dmnd_db) from diamond_db
	output:
		file('diamond_results.tsv')

	script:
	"""
	diamond blastx --threads $threads --db $dmnd_db --query $assembly_fasta --out diamond_results.tsv --outfmt 6 qseqid qlen qstart qend sseqid slen sstart send evalue score length pident nident mismatch gaps qframe cigar sscinames --header --index-chunks 1
	"""
}


def help() {
    println ""
    println "Program: denovo.nf"
    println "Developer: Steven Lakin, USDA APHIS"
    println "Description: De novo assembly and BLAST pipeline for short paired-end, hybrid short/long, or Nanopore long data."
    println ""
    println "Usage: denovo.nf [options]"
    println "Example: denovo.nf"
    println ""
    println "Input/output options:"
    println ""
    println "    --forward       STR      path to forward FASTQ input or concatenated Nanopore FASTQ input for a single sample"
    println "    --reverse       STR      path to reverse FASTQ input for a single sample"
    println "    --hybrid        STR      optional path to concatenated Nanopore FASTQ input for hybrid assembly"
    println "    --output        STR      path to output directory"
    println "    --host          STR      optional path to host genome FASTA file for host read depletion"
    println "    --reference     STR      optional path to target genome FASTA file for calculating assembly metrics"
    println "    --diamond_db    STR      path to indexed Diamond tBLASTx database"
    println "    --scratch_dir   STR      optional path to large working directory for scratch files"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of process threads, default 1 (max thread use = maxForks * threads)"
    println "    --metagenomic   FLAG     use metaSPAdes or Flye for metagenomic sample (cannot be used with hybrid assembly)"
    println "    --nanopore      FLAG     the input reads are concatenated Nanopore FASTQ data only (uses Flye assembler)"
    println "    --use_index     FLAG     use the pre-built BWA indices for host reference genome"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}
