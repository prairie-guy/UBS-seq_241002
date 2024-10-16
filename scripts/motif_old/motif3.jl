#!/usr/bin/env julia
using BioSequences
using FASTX
using CodecZlib

function extract_motifs(input_file, output_file, motif_length, reference_fasta)
    # Load the reference genome
    genome = Dict{String, LongDNA{4}}()
    open(FASTA.Reader, reference_fasta) do reader
        for record in reader
            genome[FASTA.identifier(record)] = LongDNA{4}(FASTA.sequence(LongDNA{4}, record))
        end
    end
    # Process the gzipped file directly using CodecZlib
    open(output_file, "w") do out
        open(GzipDecompressorStream, input_file) do file
            # Read and write the header, adding the new "Motif" column
            header = readline(file)
            println(out, "$header\tMotif")

            # Process each line
            for line in eachline(file)
                fields = split(line, '\t')
                if length(fields) != 9
                    continue
                end
                sample, chrom, pos, strand = fields[1], fields[2], parse(Int, fields[3]), fields[4]
                
                # Extract the sequence
                seq = try
                    chr_seq = genome[chrom]
                    if strand == "+"
                        if pos < 1 || pos + motif_length - 1 > length(chr_seq)
                            throw(BoundsError(chr_seq, pos:pos+motif_length-1))
                        end
                        chr_seq[pos:pos+motif_length-1]
                    else
                        if pos - motif_length + 1 < 1 || pos > length(chr_seq)
                            throw(BoundsError(chr_seq, pos-motif_length+1:pos))
                        end
                        reverse_complement(chr_seq[pos-motif_length+1:pos])
                    end
                catch e
                    dna"N"^motif_length
                end

                # Write output to file
                println(out, "$line\t$seq")
            end
        end
    end
end
# Main execution
if length(ARGS) != 4
    println(stderr, "Usage: julia script.jl <input.tsv.gz> <output_file> <motif_length> <reference_fasta>")
    exit(1)
end

input_file = ARGS[1]
output_file = ARGS[2]
motif_length = parse(Int, ARGS[3])
reference_fasta = ARGS[4]

extract_motifs(input_file, output_file, motif_length, reference_fasta)
