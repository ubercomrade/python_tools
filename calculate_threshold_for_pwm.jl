import CSV
import DataFrames
import Random
using ArgParse
using Distributed


function read_fasta(path)
    container = String[]
    for line in eachline(path)
        if '>' != line[1]
            line = uppercase(line)
            container = push!(container, line)
            rc_line = reverse_complement(line)
            container = push!(container, rc_line)
        else
            continue
        end
    end
    return(container)
end


@everywhere function calculate_score(site::String, pwm::Dict{Char,Array{Float64, 1}})
    score = 0.0
    len = length(site)
    for (index, nuc) in enumerate(site)
        score += pwm[site[index]][index]
    end
    return(score::Float64)
end


function reverse_complement(site::String)
    complement = Char[]
    for i in site
        if i == 'C'
            complement = push!(complement, 'G')
        elseif i == 'G'
            complement = push!(complement, 'C')
        elseif i == 'A'
            complement = push!(complement, 'T')
        elseif i == 'T'
            complement = push!(complement, 'A')
        elseif i == 'N'
            complement = push!(complement, 'N')
        end
    end
    return(join(reverse(complement)))
end


function read_pwm(path::String)
    pwm = Dict('A' => Float64[],
    'C' => Float64[],
    'G' => Float64[],
    'T' => Float64[])

    open(path) do file
        for ln in eachline(file)
            if ln[1] == '>'
                continue
            else
                ln = split(ln)
                pwm['A'] = push!(pwm['A'], parse(Float64, ln[1]))
                pwm['C'] = push!(pwm['C'], parse(Float64, ln[2]))
                pwm['G'] = push!(pwm['G'], parse(Float64, ln[3]))
                pwm['T'] = push!(pwm['T'], parse(Float64, ln[4]))
            end
        end
    end
    return(pwm)
end


function scan_peak(peak::String, len_of_site::Int64, pwm::Dict{Char,Array{Float64, 1}})
    scores = Float64[]
    len = length(peak)
    for i in 1:len - len_of_site
      site = peak[i:i + len_of_site - 1]
      if 'N' in site
          continue
      end
      scores = push!(scores, calculate_score(site, pwm))
    end
    return(scores)
end

function scan_peaks(peaks::Array{String, 1}, len_of_site::Int64, pwm::Dict{Char,Array{Float64, 1}})
    l = length(peaks)
    scores = Array{Array{Float64, 1}, 1}()
    Threads.@threads for peak in peaks
        scores =  push!(scores, scan_peak(peak::String, len_of_site::Int64, pwm::Dict{Char,Array{Float64, 1}}))
    end
    return(reduce(vcat, scores::Array{Array{Float64, 1}, 1}))
end


function calculate_thresholds(peaks::Array{String, 1}, pwm::Dict{Char,Array{Float64, 1}}, len_of_site::Int64)

    # scores = broadcast(peak::String -> scan_peak(peak::String, len_of_site::Int64, pwm::Dict{Char,Array{Float64, 1}}), peaks)
    # scores = reduce(vcat, scores::Array{Array{Float64, 1}, 1});

    scores = scan_peaks(peaks, len_of_site::Int64, pwm::Dict{Char,Array{Float64, 1}});
    scores = sort!(scores::Array{Float64, 1}, rev=true)

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
        "pwm"
            help = "path to .pwm file"
            required = true
        "output"
            help = "path to write results"
            required = true

    end
    return parse_args(s)
end


function main()
    args = parse_commandline()
    fasta_path = args["fasta"]
    pwm_path = args["pwm"]
    output = args["output"]

    pwm = read_pwm(pwm_path);
    len_of_site = length(pwm['A']);
    peaks = read_fasta(fasta_path);

    res = calculate_thresholds(peaks, pwm, len_of_site)
    CSV.write(output, res, delim='\t')

end

#main()

# pwm_path = "/Users/anton/Google Диск/PhD/Расчеты/CHOSEN-TFS-SCAN-5000/AR_41593/MOTIFS/AR_41593_OPTIMAL_MOTIF.pwm"
# fasta_path = "/Users/anton/Documents/DATA/PROMOTERS/mm10_promoters.fa"
#
# pwm = read_pwm(pwm_path);
# len_of_site = length(pwm['A']);
# peaks = read_fasta(fasta_path);
# @time s = scan_peaks(peaks, len_of_site::Int64, pwm::Dict{Char,Array{Float64, 1}});
