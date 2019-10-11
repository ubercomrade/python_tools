import CSV
import DataFrames
import Random
using ArgParse


function read_logOddsZoops(path)
    container = Array{Char, 1}[]
    index = 0
    for line in eachline(path)
        if index == 0
            index += 1
            continue
        end
        container = push!(container, collect(split(line)[5]))
    end
    container = check_n(container)
    return(container)
end


function check_n(sites::Array{Array{Char, 1},1})
    indexes = Bool[]
    for site in sites
        indexes = push!(indexes, !('N' in site))
    end
    return(sites[indexes])
end


function read_sites(path)
    container = Array{Char, 1}[]
    for line in eachline(path)
        if '>' != line[1]
            container = push!(container, collect(line))
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


function calculate_score(site::Array{Char, 1}, bamm::Dict{String,Array{Float64, 1}}, order::Int64)
    score = 0.0
    l = length(site)
    for index in 1:order
        score += bamm[join(site[1:index])][index]
    end
    for index in 1:l - order
        score += bamm[join(site[index:index+order])][index + order]
    end
    return(score)
end


function calculate_scores(sites::Array{Array{Char, 1}}, bamm::Dict{String,Array{Float64, 1}}, order::Int64)
    scores = Float64[]
    for site in sites
        score = calculate_score(site, bamm, order)
        scores = push!(scores, score)
    end
    return(scores)
end


function calculate_bamm_model(dir::String, order::Int64)
    fasta_path = string(dir, "/train.fasta")
    sites_path = string(dir, "/sites.txt")
    #while ! success(`BaMMmotif $dir $fasta_path --bindingSiteFile $sites_path --EM`)
    #end
    run(`BaMMmotif $dir $fasta_path --bindingSiteFile $sites_path --EM`)
    bamm_path = string(dir, "/train_motif_1.ihbcp")
    bg_path = string(dir, "/train.hbcp")
    bamm = read_bamm(bamm_path, bg_path, order)
    return(bamm)
end


function create_train_fasta_and_sites(sites::Array{Array{Char, 1}}, dir::String)
    open(string(dir,"/train.fasta"), "w") do file
        for (index, site) in enumerate(sites)
            write(file, string(">",index,'\n',join(site),'\n'))
        end
    end

    open(string(dir,"/sites.txt"), "w") do file
        for site in sites
            write(file, string(join(site),'\n'))
        end
    end
end


function shuffling(sites::Array{Array{Char, 1}})
    container = Array{Array{Char, 1},1}()
    s = join(join.(sites));
    r = Random.shuffle(collect(s));
    l = length(r)
    step = length(sites[1])
    for i in 1:step:l
        container = push!(container, r[i:i+step-1])
    end
    return(container)
end


function bootstrap_bamm(sites::Array{Array{Char, 1}}, size_of::Int64, order::Int64, dir::String)
    true_scores = Float64[]
    false_scores = Float64[]
    number_of_sites = length(sites)

    for i in 1:10
        index_train = Random.randsubseq(1:number_of_sites, 0.9)
        index_test = setdiff(1:number_of_sites, index_train)
        index_shuffle = Random.rand(1:number_of_sites, size_of)

        create_train_fasta_and_sites(sites[index_train], dir)
        bamm = calculate_bamm_model(dir, order)
        true_scores = vcat(true_scores, calculate_scores(sites[index_test], bamm, order))
        false_scores = vcat(false_scores, calculate_scores(
                Random.shuffle.(sites[index_shuffle]),
                bamm, order))
    end

    true_scores = sort(true_scores, rev=true)
    false_scores = sort(false_scores, rev=true)

    tpr = Float64[]
    tpr_actual = Float64[]
    fpr = Float64[]
    scores = Float64[]

    false_length = length(false_scores)
    true_length = length(true_scores)

    for i in 0.05:0.05:1.0
        s = true_scores[Int(round(true_length * i))]
        tpr_actual = push!(tpr_actual, sum(true_scores .>= s) / true_length)
        tpr = push!(tpr, i)
        fpr = push!(fpr, sum(false_scores .>= s) / false_length)
        scores = push!(scores, s)
    end

    df = DataFrames.DataFrame(Scores = scores, TPR = tpr, ACTUAL_TPR = tpr_actual, FPR = fpr)
    return(df)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "output"
            help = "path to write results"
            required = true
        "input"
            help = "path to .logOddsZoops file from BaMMmotif output with flags --scoreSeqset and --saveLogOdds"
            required = true
        "--size", "-s"
            help = "size of negative sites on each bootstrap iteration"
            arg_type = Int
            default = 100000
        "--order", "-o"
            help = "model order"
            required = false
            arg_type = Int
            default = 2
        "--tmp", "-t"
            help = "dir for tmp files"
            required = false
            arg_type = String
            default = "./tmp"

    end
    return parse_args(s)
end


function main()
    args = parse_commandline()
    path = args["input"]
    out = args["output"]
    size_of = args["size"]
    order = args["order"]
    tmp_dir = args["tmp"]

    if !isdir(tmp_dir)
        mkdir(tmp_dir)
    end

    sites = read_logOddsZoops(path)
    df = bootstrap_bamm(sites, size_of, order, tmp_dir)
    CSV.write(out, df, delim='\t')

    rm(tmp_dir, recursive=true)
end

main()
