#!/usr/bin/env julia

using BioSequences
using FASTX

function extract_motifs(input_file, output_file, motif_length, reference_fasta; lines=10000000)
    # Load the reference genome
    genome = Dict{String, LongDNA{4}}()
    open(FASTA.Reader, reference_fasta) do reader
        for record in reader
            genome[FASTA.identifier(record)] = LongDNA{4}(FASTA.sequence(LongDNA{4}, record))
        end
    end

    # Decompress the input file using pigz
    decompressed_file = tempname()
    run(pipeline(`pigz -d -k -c $input_file`, stdout=decompressed_file))

    # Process the decompressed file and write to output file
    open(output_file, "w") do out
        open(decompressed_file) do file
            # Read and write the header, adding the new "Motif" column
            header = readline(file)
            println(out, "$header\tMotif")

            # Initialize buffer
            buffer = IOBuffer()
            line_count = 0

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

                # Write output to buffer
                println(buffer, "$line\t$seq")
                line_count += 1

                # If buffer is full, write to file and reset
                if line_count >= lines
                    write(out, take!(buffer))
                    line_count = 0
                end
            end

            # Write any remaining lines in the buffer
            if line_count > 0
                write(out, take!(buffer))
            end
        end
    end

    # Clean up the decompressed file
    rm(decompressed_file)
end

# Main execution
if length(ARGS) < 4 || length(ARGS) > 5
    println(stderr, "Usage: julia script.jl <input.tsv.gz> <output_file> <motif_length> <reference_fasta> [buffer_lines]")
    exit(1)
end

input_file = ARGS[1]
output_file = ARGS[2]
motif_length = parse(Int, ARGS[3])
reference_fasta = ARGS[4]
buffer_lines = length(ARGS) == 5 ? parse(Int, ARGS[5]) : 10000000

extract_motifs(input_file, output_file, motif_length, reference_fasta, lines=buffer_lines)
