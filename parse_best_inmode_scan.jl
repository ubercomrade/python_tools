import CSV
import DataFrames
import Random
using ArgParse


function read_fasta(path)
    container = Dict("ID" => String[],
    "chr" => String[],
    "start" => Int64[],
    "end" => Int64[],
    "strand" => Char[],
    "seq" => String[])
    for line in eachline(path)
        if line[1] == '>'
            container["ID"] = push!(container["ID"], split(line, "::")[1][2:end])
            second_part = split(line, "::")[2]
            container["chr"] = push!(container["chr"], split(second_part, ":")[1])
            container["start"] = push!(container["start"], parse(Int64, split(split(second_part, ":")[2][1:end-2], "-")[1]))
            container["end"] = push!(container["end"], parse(Int64, split(split(second_part, ":")[2][1:end-2], "-")[2]))
            container["strand"] = push!(container["strand"], '+')
        else
            container["seq"] = push!(container["seq"], uppercase(line))
        end
    end
    return(container)
end

function read_best_inmode(path)
    ID = 0
    score = 0
    save_index = 0
    container = Dict("ID" => Int64[],
    "start" => Int64[],
    "end" => Int64[],
    "chr" => Char[],
    "score" => Float64[] )
    for line in eachline(path)
        line_id = parse(Int64, split(line)[1])
        line_start = parse(Int64, split(line)[2])
        line_end = parse(Int64, split(line)[3])
        line_chr = parse(Char, split(line)[4])
        line_score = parse(Float64, split(line)[5])
