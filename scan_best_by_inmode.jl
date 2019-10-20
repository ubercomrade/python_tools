import CSV
import DataFrames
import Random
using ArgParse


function read_fasta(path)
    container = Dict("ID" => String[],
    "chr" => String[],
    "start" => Int64[],
    "end" => Int64[],
    "strand" => String[],
    "seq" => String[])
    for line in eachline(path)
        if line[1] == '>'
            container["ID"] = push!(container["ID"], split(line, "::")[1][2:end])
            second_part = split(line, "::")[2]
            container["chr"] = push!(container["chr"], split(second_part, ":")[1])
            container["start"] = push!(container["start"], parse(Int64, split(split(second_part, ":")[2][1:end-2], "-")[1]))
            container["end"] = push!(container["end"], parse(Int64, split(split(second_part, ":")[2][1:end-2], "-")[2]))
            container["strand"] = push!(container["strand"], "+")
        else
            container["seq"] = push!(container["seq"], uppercase(line))
        end
    end
    return(container)
end

function read_best_inmode(path)
    ID = 0
    score = 0.0
    container = Dict("ID" => Int64[],
    "start" => Int64[],
    "end" => Int64[],
    "strand" => String[],
    "score" => Float64[] )
    
    line_id = Int64
    line_start = Int64
    line_end = Int64
    line_strand = String
    line_score = Float64

    for line in eachline(path)
        if (parse(Int64, split(line)[1]) == ID) & (parse(Float64, split(line)[5]) > score)
            line_id = parse(Int64, split(line)[1])
            line_start = parse(Int64, split(line)[2])
            line_end = parse(Int64, split(line)[3])
            line_strand = split(line)[4]
            line_score = parse(Float64, split(line)[5])
            score = line_score
            
        elseif (parse(Int64, split(line)[1]) != ID)
            
            container["ID"] = push!(container["ID"], line_id)
            container["start"] = push!(container["start"], line_start)
            container["end"] = push!(container["end"], line_end)
            container["strand"] = push!(container["strand"], line_strand)
            container["score"] = push!(container["score"], line_score)
            
            line_id = parse(Int64, split(line)[1])
            ID = line_id
            line_start = parse(Int64, split(line)[2])
            line_end = parse(Int64, split(line)[3])
            line_strand = split(line)[4]
            line_score = parse(Float64, split(line)[5])
            score = line_score
        end
    end
    return(container)
end


function reverse_complement(site::String)
    complement = ""
    for i in site
        #println(i)
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


function create_bed(fasta, inmode)
    container = Dict("ID" => String[],
    "chr" => String[],
    "start" => Int64[],
    "end" => Int64[],
    "strand" => String[],
    "site" => String[],
    "score" => Float64[])

    for (score, start, enD, strand, ID) in zip(inmode["score"], inmode["start"], inmode["end"], inmode["strand"], inmode["ID"])
        if strand == "-"
            container["site"] = push!(container["site"], reverse_complement(fasta["seq"][ID+1][start+1:enD]))
        else
            container["site"] = push!(container["site"], fasta["seq"][ID+1][start+1:enD])
        end

        container["chr"] = push!(container["chr"], fasta["chr"][ID+1])
        container["start"] = push!(container["start"], fasta["start"][ID+1] + start)
        container["end"] = push!(container["end"], fasta["start"][ID+1] + enD)   
        container["ID"] = push!(container["ID"], string("peaks_", string(ID)))
        container["strand"] = push!(container["strand"], strand)
        container["score"] = push!(container["score"], score)
    end
    
    df = DataFrames.DataFrame(chr = container["chr"],
        st = container["start"],
        en = container["end"],
        id = container["ID"],
        score = container["score"],
        strand = container["strand"],
        site = container["site"])
    
    return(df)
end


function scan_all_by_inmode(path_to_inmode, path_to_model, path_to_fasta, tmp_dir)
    #run(`/home/anton/Programs/jdk-9/bin/java -Xmx4096m -Xms1024m --add-modules java.xml.bind -jar $path_to_inmode scan
    #run(`java -jar $path_to_inmode scan
    run(`/Users/anton/Documents/Programs/jre-9.0.4.jre/Contents/Home/bin/java -Xmx4096m -Xms1024m --add-modules java.xml.bind -jar $path_to_inmode scan
    i="$path_to_model"
    id="$path_to_fasta"
    f=1.0
    outdir=$tmp_dir
    bs=true`);
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "output"
            help = "path to write results"
            required = true
        "input"
            help = "path to xml file InMoDe model"
            required = true
        "fasta"
            help = "path to fasta file with sites"
            required = true
        "inmode"
            help = "path to inmode program"
            required = true
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
    path_to_model = args["input"]
    out = args["output"]
    path_to_fasta = args["fasta"]
    tmp_dir = args["tmp"]
    path_to_inmode = args["inmode"]
    
    if !isdir(tmp_dir)
        mkdir(tmp_dir)
    end

    scan_all_by_inmode(path_to_inmode, path_to_model, path_to_fasta, tmp_dir)
    fasta = read_fasta(path_to_fasta)
    inmode_path = string(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED")
    inmode = read_best_inmode(inmode_path)
    df = create_bed(fasta, inmode)
    rm(tmp_dir, recursive=true)
    CSV.write(out, df, delim='\t',writeheader=false)
end

main()