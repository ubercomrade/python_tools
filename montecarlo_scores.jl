import Random
using ArgParse
using Statistics
using StatsPlots
# using Plots

function read_list(path)
    container = Float64[]
    open(path, "r") do file
        for ln in eachline(file)
            container = push!(container, parse(Float64,ln))
        end
    end
    return(container)
end


function sort_scores(scores, lengths)
    model_tuple = Tuple[]
    for (s, l) in zip(scores, lengths)
        model_tuple = push!(model_tuple, (s, l))
    end

    score_sorted = Float64[]
    model_tuple = sort(model_tuple, by = tuple -> last(tuple))
    for i in sort(model_tuple, by = tuple -> last(tuple))
        score_sorted = push!(score_sorted, i[1])
    end
    return(score_sorted)
end


function create_bins(len,step,shift)

    maxLen = maximum(len);
    maxLen = maxLen + 10 - maxLen % 10;
    minLen = minimum(len);
    minLen = minLen - minLen % 10;

    bins_1 = Array{Array{Int64, 1}, 1}()
    flag = false
    prev = 0
    total_sum = 0
    for i in minLen:step:maxLen
        s = sum((len .< i + step) .& (len .>= i - prev))
        if s >= 10
            bins_1 = push!(bins_1, collect(1+total_sum:s+total_sum))
            prev = 0
            total_sum += s
        else
            prev += step
        end
    end

    bins_2 = Array{Array{Int64, 1}, 1}()
    flag = false
    prev = 0
    bins_2 = push!(bins_2, collect(1:sum(len .< minLen+shift)))
    total_sum = sum(len .< minLen+shift)
    for i in minLen+shift:step:maxLen
        s = sum((len .< i + step) .& (len .>= i - prev))
        if s >= 10
            bins_2 = push!(bins_2, collect(1+total_sum:s+total_sum))
            prev = 0
            total_sum += s
        else
            prev += step
        end
    end

    bins_2[end] = filter(t -> t <= bins_1[end][end], bins_2[end])

    return(bins_1, bins_2)
end


function creat_distribution_by_montecarlo(first_model_scores::Array{Float64, 1}, second_model_scores::Array{Float64, 1},
        first_model_threshold::Float64, second_model_threshold::Float64,
        iterations::Int64, bins_1::Array{Array{Int64,1},1}, bins_2::Array{Array{Int64,1},1})
    container = Float64[]
    for i in 1:iterations
        first_model_scores_shuffled_1 = Float64[]
        second_model_scores_shuffled_1 = Float64[]

        for i in bins_1
            first_model_scores_shuffled_1 = vcat(first_model_scores_shuffled_1, first_model_scores[Random.shuffle(i)])
            second_model_scores_shuffled_1 = vcat(second_model_scores_shuffled_1, second_model_scores[Random.shuffle(i)])
        end

        first_model_scores_shuffled_2 = Float64[]
        second_model_scores_shuffled_2 = Float64[]

        for i in bins_2
            first_model_scores_shuffled_2 = vcat(first_model_scores_shuffled_2, first_model_scores_shuffled_1[Random.shuffle(i)])
            second_model_scores_shuffled_2 = vcat(second_model_scores_shuffled_2, second_model_scores_shuffled_1[Random.shuffle(i)])
        end

        container = push!(container, sum((first_model_scores_shuffled_2 .> first_model_threshold) .& (second_model_scores_shuffled_2 .> second_model_threshold)))
    end
    return(container)
end


function plot_scores(first_model_scores, second_model_scores,
        first_model_threshold, second_model_threshold, path_out_plot)
    scatter(first_model_scores, second_model_scores, leg=false, xlabel="PWM", ylabel="BAMM")
    vline!([first_model_threshold], linewidth=2, linecolor=:red)
    hline!([second_model_threshold], linewidth=2, linecolor=:red)
    savefig(path_out_plot)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "firstModel"
            help = "path to first model scores"
            required = true
        "secondModel"
            help = "path to second model scores"
            required = true
        "firstThr"
            help = "threshold for first model"
            required = true
            arg_type = Float64
        "secondThr"
            help = "threshold for second model"
            required = true
            arg_type = Float64
        "length"
            help = "path to file with length of peaks"
            required = true
        "--out", "-o"
            help = "path to write scatter plot"
            required = false
    end
    return parse_args(s)
end


function main()
    args = parse_commandline()
    first_model_path = args["firstModel"]
    second_model_path = args["secondModel"]
    length_of_peaks_path = args["length"]
    first_model_threshold = args["firstThr"]
    second_model_threshold = args["secondThr"]
    path_out_plot = args["out"]


    first_model_scores = read_list(first_model_path);
    second_model_scores = read_list(second_model_path);
    len = read_list(length_of_peaks_path);

    step = 50;
    shift = 25;
    iterations = 10000;
    bins_1, bins_2 = create_bins(len,step,shift);

    first_model_scores = sort_scores(first_model_scores, len);
    second_model_scores = sort_scores(second_model_scores, len);
    len = sort(len);

    common_scores = sum((first_model_scores .> first_model_threshold) .& (second_model_scores .> second_model_threshold));
    container = creat_distribution_by_montecarlo(first_model_scores,
        second_model_scores,
        first_model_threshold,
        second_model_threshold,
        iterations, bins_1, bins_2);

    rand = mean(container)
    sd = std(container)
    zscore = (common_scores - mean(container)) / std(container)

    if path_out_plot != nothing
        plot_scores(first_model_scores, second_model_scores,
            first_model_threshold, second_model_threshold, path_out_plot)
        println("$common_scores\t$rand\t$sd\t$zscore\n")
    else
        println("$common_scores\t$rand\t$sd\t$zscore\n")
    end
end

main()
