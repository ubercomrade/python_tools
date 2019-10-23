import CSV
import DataFrames
import Random
using ArgParse
using Distributed


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


function calculate_score(site::String, pwm::Dict{Char,Array{Float64, 1}})
    score = 0.0
    for (index, nuc) in enumerate(site)
        score += pwm[nuc][index]
    end
    return(score)
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


function scan_sites(sites::Array{String, 1}, pwm)
    scores = Float64[]
    for i in sites
        scores = push!(scores, calculate_score(i, pwm))
    end
    return(scores)
end


function split_to_sites_atom(peak::String, len_of_site::Int)
    sites = String[]
    for i in 1:length(peak) - len_of_site
        site = peak[i:i + len_of_site - 1]
        if 'N' in site
            continue
        end
        sites = push!(sites, site)
        sites = push!(sites, reverse_complement(site))
    end
    return(sites)
end


function calculate_thresholds(peaks::Array{String, 1}, pwm::Dict{Char,Array{Float64, 1}}, len_of_site::Int64)
    scores = Float64[]
    
#    for peak in peaks
#        for i in 1:length(peak) - len_of_site
#            site = peak[i:i + len_of_site - 1]
#            if 'N' in site
#                continue
#            end

#            score = calculate_score(site, pwm)
#            scores = push!(scores, score)

#            score = calculate_score(reverse_complement(site), pwm)
#            scores = push!(scores, score)

#        end
#    end
    
    for peak in peaks
          sites = String[]
          for i in 1:length(peak) - len_of_site
              site = peak[i:i + len_of_site - 1]
              if 'N' in site
                  continue
              end
              sites = push!(sites, site)
              sites = push!(sites, reverse_complement(site))
          end

          for i in sites
              scores = push!(scores, calculate_score(i, pwm))
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
    fasta = read_fasta(fasta_path);
    len_of_site = length(pwm['A']);

    res = calculate_thresholds(fasta, pwm, len_of_site)
    CSV.write(output, res, delim='\t')
    
end

main()
