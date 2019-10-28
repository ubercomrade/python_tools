import CSV
import DataFrames
import Random
using ArgParse

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


function create_pcm(sites)
    l = length(sites[1])
    pcm = Dict{Char, Array{Float64, 1}}(
        'A' => zeros(Float64, l),
        'C' => zeros(Float64, l),
        'G' => zeros(Float64, l),
        'T' => zeros(Float64, l))

    for i in sites, (index, nuc) in enumerate(i)
        pcm[nuc][index] += 1.0
    end
    return(pcm)
end


function create_pfm(pcm::Dict{Char,Array{Float64, 1}}, number_of_sites::Int64)
    pfm = Dict{Char, Array{Float64, 1}}()
    for nuc in keys(pcm)
        pfm[nuc] = (pcm[nuc] .+ 1.0/number_of_sites) ./ (number_of_sites + 1.0)
    end
    return(pfm)
end


function create_pwm(pfm::Dict{Char,Array{Float64, 1}}, background::Dict{Char,Float64})
    pwm = Dict{Char, Array{Float64, 1}}()
    for nuc in keys(pfm)
        pwm[nuc] = log2.(pfm[nuc] ./ background[nuc])
    end
    return(pwm)
end


function calculate_background(sites)
    all_nucleotides = join(sites)
    l = length(all_nucleotides)
    background = Dict{Char, Float64}()
    for nuc in all_nucleotides
        if !(nuc in keys(background))
            background[nuc] = 1
        else
            background[nuc] += 1
        end
    end

    for nuc in keys(background)
        background[nuc] = background[nuc] / l
    end
    return(background)
end


function make_pwm(sites)
    number_of_sites = length(sites)
    pcm = create_pcm(sites)
    pfm = create_pfm(pcm, number_of_sites)
    background = calculate_background(sites)
    pwm = create_pwm(pfm, background)
    return(pwm::Dict{Char,Array{Float64, 1}})
end


function calculate_score(site, pwm::Dict{Char,Array{Float64, 1}})
    score = 0.0
    for (index, nuc) in enumerate(site)
        score += pwm[nuc][index]
    end
    return(score)
end


function calculate_scores(sites, pwm::Dict{Char,Array{Float64, 1}})
    scores = Float64[]
    for site in sites
        score = calculate_score(site, pwm)
        scores = push!(scores, score)
    end
    return(scores)
end


function bootstrap_pwm(sites, size_of::Int64)
    true_scores = Float64[]
    false_scores = Float64[]
    number_of_sites = length(sites)

    for i in 1:10
        index_train = Random.randsubseq(1:number_of_sites, 0.9)
        index_test = setdiff(1:number_of_sites, index_train)
        index_shuffle = Random.rand(1:number_of_sites, size_of)

        pwm = make_pwm(sites[index_train])
        true_scores = vcat(true_scores, calculate_scores(sites[index_test], pwm))
        false_scores = vcat(false_scores, calculate_scores(
                Random.shuffle.(sites[index_shuffle]),
                pwm))
    end

    println(typeof(false_scores[1]))

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
            help = "path to sites of PWM model"
            required = true
        "--size", "-s"
            help = "size of negative sites on each bootstrap iteration"
            arg_type = Int
            default = 100000

    end
    return parse_args(s)
end


function main()
    args = parse_commandline()
    path = args["input"]
    out = args["output"]
    size_of = args["size"]
    sites = read_sites(path)
    df = bootstrap_pwm(sites, size_of)
    CSV.write(out, df, delim='\t')
end

# function main()
#     s = parse_commandline()
#     println("Parsed args:")
#     println(s)
#     println(s["input"])
#     println(s["output"])
#     println(s["size"])
#     for (arg,val) in s
#         println("  $arg  =>  $val")
#     end
# end

main()
