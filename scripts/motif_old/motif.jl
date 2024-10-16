#!/usr/bin/env julia

using BioSequences
using FASTX
using CodecZlib

function extract_motifs(input_file, motif_length, reference_fasta)
    # Load the reference genome
    genome = Dict{String, LongDNA{4}}()
    reader = open(FASTA.Reader, reference_fasta)
    for record in reader
        genome[FASTA.identifier(record)] = LongDNA{4}(FASTA.sequence(LongDNA{4}, record))
    end
    close(reader)

    # Open the gzipped input file
    open(GzipDecompressorStream, input_file) do file
        # Read and print the header, adding the new "Motif" column
        header = readline(file)
        println("$header\tMotif")

        # Process each line
        for line in eachline(file)
            fields = split(line, '\t')
            if length(fields) != 9
                println(stderr, "Warning: Unexpected number of fields in line: $line")
                continue
            end

            sample, chrom, pos, strand = fields[1], fields[2], parse(Int, fields[3]), fields[4]

            # Extract the sequence
            try
                chr_seq = genome[chrom]
                if pos < 1 || pos + motif_length - 1 > length(chr_seq)
                    throw(BoundsError(chr_seq, pos:pos+motif_length-1))
                end
                seq = if strand == "+"
                    chr_seq[pos:pos+motif_length-1]
                else
                    reverse_complement(chr_seq[pos-motif_length+1:pos])
                end
                println("$line\t$seq")
            catch e
                if isa(e, KeyError)
                    println(stderr, "Warning: Chromosome $chrom not found in reference genome.")
                elseif isa(e, BoundsError)
                    println(stderr, "Warning: Position out of bounds for chromosome $chrom (position: $pos, chromosome length: $(length(get(genome, chrom, dna"")))).")
                else
                    println(stderr, "Unexpected error: $e")
                end
                println("$line\t$(dna"N"^motif_length)")
            end
        end
    end
end

# Main execution
if length(ARGS) != 3
    println("Usage: julia script.jl <input.tsv.gz> <motif_length> <reference.fasta>")
    exit(1)
end

input_file = ARGS[1]
motif_length = parse(Int, ARGS[2])
reference_fasta = ARGS[3]

extract_motifs(input_file, motif_length, reference_fasta)
