

module MergeSampleFname
using CSV
using DataFrames
using Glob
export merge_sample_fname
"""
    merge_sample_fname(in_csv::String, dir_path::String, out_csv::String="sample_fname.csv")

Merge sample information from a CSV file with filename data from a specified directory.

# Arguments
- `in_csv::String`: Path to the input CSV file containing sample information.
                    This file must include 'SampleID' and 'Content' columns.
- `dir_path::String`: Path to the directory containing the .fastq.gz or .fq.gz files.
                      Sequences in this directory are expected to contain the SampleID.
- `out_csv::String`: Path where the output CSV file will be written.
                     Defaults to 'sample_fname.csv' in the current directory.

# Returns
Nothing. The function writes the results to a CSV file and prints a confirmation message.
"""
function merge_sample_fname(in_csv::String, dir_path::String, out_csv::String="sample_fname.csv")
    df = CSV.read(in_csv, DataFrame)
    files = glob("*.[ff][aq]*.gz", dir_path)
    files_data = DataFrame(SampleID = String[], Std = String[], Sequence = String[])
    for file in files
        m = match(r"-(.+?)_.*_(R[12])_", basename(file))
        if !isnothing(m)
            push!(files_data, (m.captures[1], m.captures[2], string(file)))
        end
    end
    result = innerjoin(df, files_data, on=:SampleID)
    result.FullID = result.SampleID .* "_" .* result.Std
    result = result[:, [:SampleID, :Std, :FullID, :Content, :Sequence]]
    CSV.write(out_csv, result)
    println("Results written to $out_csv")
end
end
