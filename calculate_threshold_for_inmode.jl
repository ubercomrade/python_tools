import CSV
import Random
using ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "fasta"
            help = "path to fasta file"
            required = true
        "input"
            help = "path to xml file (InMoDe model)"
            required = true
        "inmode"
            help = "path to InMoDe java program"
            required = true
        "output"
            help = "path to write results"
            required = true
        "--tmp", "-t"
            help = "dir for tmp files"
            required = false
            arg_type = String
            default = "./tmp"
        "--java", "-j"
            help = "path to java"
            required = false
            arg_type = String
            default = "java"


    end
    return parse_args(s)
end


function calculate_scores(path_to_inmode::String, path_to_model::String,
    path_to_fasta::String, path_to_java::String, tmp_dir::String)
    container = Float64[]
    run(`$path_to_java -Xmx4096m -Xms1024m --add-modules java.xml.bind -jar $path_to_inmode scan i="$path_to_model" id="$path_to_fasta" f=1.0 outdir=$tmp_dir bs=true`);
    for line in eachline(string(tmp_dir,"/Motif_hits_from_SequenceScan(1.0).BED"))
        container = push!(container, parse(Float64, split(line)[5]))
    end
    run(`rm $tmp_dir/Binding_sites_from_SequenceScan\(1.0\).txt $tmp_dir/Motif_hits_from_SequenceScan\(1.0\).BED $tmp_dir/protocol_scan.txt $tmp_dir/$tag.fa`);
    return(container)
end



function calculate_thresholds(path_to_inmode::String, path_to_model::String,
    path_to_fasta::String, path_to_java::String, tmp_dir::String)

    scores = calculate_scores(path_to_inmode, path_to_model, path_to_fasta, path_to_java, tmp_dir)
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


function main()
    args = parse_commandline()
    path_to_model = args["input"]
    out = args["output"]
    path_to_fasta = args["fasta"]
    tmp_dir = args["tmp"]
    path_to_inmode = args["inmode"]
    path_to_java = args["java"]

    if !isdir(tmp_dir)
        mkdir(tmp_dir)
    end

    res = calculate_thresholds(path_to_inmode, path_to_model, path_to_fasta, path_to_java, tmp_dir)
    CSV.write(output, res, delim='\t')

    rm(tmp_dir, recursive=true)

end

main()
