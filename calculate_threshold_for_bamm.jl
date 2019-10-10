import CSV
import DataFrames
import Random
using ArgParse


function read_fasta(path)
    container = String[]
    for line in eachline(path)
        if '>' != line[1]
            container = push!(container, (uppercase(line)))
        else
            continue
        end
    end
    return(container)
end


function product_of_nucleotides(order::Int64)
    container = String[]
    if order == 0
        return(["A","C","G","T"])
    elseif order == 1
        for i in ["A","C","G","T"]
            for j in ["A","C","G","T"]
                container = push!(container, string(i, j))
            end
        end
    else
        order -= 1
        for i in ["A","C","G","T"]
            for j in product_of_nucleotides(order)
                container = push!(container, string(i, j))
            end
        end
    end
    return(container)
end


function support_read_bamm_logodds(bamm_path::String, bg_path::String, order::Int64)
    
    bamm = Dict{String, Array{Float64, 1}}()
    nucleotides = product_of_nucleotides(order)
    for i in nucleotides
        bamm[i] = Float64[]
    end
    index = 0
    for line in eachline(bamm_path)
        index += 1
        if line == ""
            index = 0
            print(line)
            continue
        elseif index == order + 1
            values = parse.(Float64, split(line))
            for (nuc, val) in zip(nucleotides, values)
                bamm[nuc] = push!(bamm[nuc], val)
            end
        end
    end
    
    index = 0
    
    bg = Dict{String, Array{Float64, 1}}()
    for i in nucleotides
        bg[i] = Float64[]
    end
    for line in eachline(bg_path)
        if line[1] == '#'
            continue
        end
        index += 1
        if index == order + 1
            values = parse.(Float64, split(line))
            for (nuc, val) in zip(nucleotides, values)
                bg[nuc] = push!(bg[nuc], val)
            end
        end
    end
    
    for nuc in nucleotides
    bamm[nuc] = log2.(bamm[nuc] ./ bg[nuc])
    end
    
    return(bamm)
end


function read_bamm(bamm_path::String, bg_path::String, order::Int64)
    bamm = support_read_bamm_logodds(bamm_path, bg_path, 0)
    for i in 1:order
        bamm = merge!(bamm, support_read_bamm_logodds(bamm_path, bg_path, i))
    end
    return(bamm)
end


function reverse_complement(site::String)
    complement = ""
    for i in site
        if i == 'C'
            complement = string('G', complement)
        elseif i == 'G'
            complement = string('C', complement)
        elseif i == 'A'
            complement = string('T', complement)
        elseif i == 'T'
            complement = string('A', complement)
        end
    end
    return(complement)
end


function calculate_score(site::String, bamm::Dict{String,Array{Float64, 1}}, order::Int64, len_of_site::Int64)
    score = 0.0
    for index in 1:order
        score += bamm[site[1:index]][index] 
    end
    for index in 1:len_of_site - order
        score += bamm[site[index:index+order]][index + order]
    end
    return(score)
end


function insert_sort(v::Vector, x)
    return((splice!(v, searchsorted(v,x), [x]); v))
end


function calculate_scores(peaks::Array{String, 1}, bamm::Dict{String,Array{Float64, 1}}, order::Int64, len_of_site::Int64)
    scores = Float64[]
    for peak in peaks
        for i in 1:length(peak) - len_of_site
            site = peak[i:i + len_of_site]
            if 'N' in site
                continue
            end
            score = calculate_score(site, bamm, order, len_of_site)
            scores = push!(scores, score)
            #scores = insert_sort(scores, score)
            score = calculate_score(reverse_complement(site), bamm, order, len_of_site)
            scores = push!(scores, score)
            #scores = insert_sort(scores, score)
        end
    end
    
    scores = sort(scores, rev=true)
    
    fpr_actual = Float64[]
    fpr = Float64[]
    scores_to_write = Float64[]

    scores_length = length(scores)

    for i in [5e-3,1e-3,5e-4,1e-4,5e-5,1e-5]
        s = scores[Int(round(scores_length * i))]
        #fpr_actual = push!(fpr_actual, sum(res .>= s) / res_length)
        fpr = push!(fpr, i)
        scores_to_write = push!(scores_to_write, s)
    end    

    df = DataFrames.DataFrame(Scores = scores_to_write, FPR = fpr)
    return(df)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "fasta"
            help = "path to fasta file"
            required = true
        "bamm"
            help = "path to .ihbcp file from BaMMmotif output"
            required = true
        "bg"
            help = "path to .hbcp file from BaMMmotif output"
            required = true
        "output"
            help = "path to write results"
            required = true
        "--order", "-o"
            help = "model order"
            required = false
            arg_type = Int
            default = 0

    end
    return parse_args(s)
end


function main()
    args = parse_commandline()
    fasta_path = args["fasta"]
    bamm_path = args["bamm"]
    bg_path = args["bg"]
    output = args["output"]
    order = args["order"]
    
    peaks = read_fasta(fasta)
    bamm = read_bamm(bamm_path, bg_path, order)
    len_of_site = length(bamm["A"])
    res = calculate_scores(peaks, bamm, order, len_of_site)
    CSV.write(output, res, delim='\t')
    
end

main()