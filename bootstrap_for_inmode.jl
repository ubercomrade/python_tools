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

function write_fasta(sites, tmp_dir, tag)
    open(string(tmp_dir, "/$tag.fa"), "w") do io
        for (index, site) in enumerate(sites)
            site = join(site)
            write(io, ">$index\n$site\n")
        end
    end
end


function calculate_scores(path_to_inmode, path_to_java::String, motif_length, tmp_dir, tag)
    container = Float64[]
    #run(`/home/anton/Programs/jdk-9/bin/java -Xmx4096m -Xms1024m --add-modules java.xml.bind -jar $path_to_inmode scan
    #run(`java -jar $path_to_inmode scan
    run(`$path_to_java -Xmx4096m -Xms1024m --add-modules java.xml.bind -jar $path_to_inmode scan
    i="$tmp_dir/Learned_DeNovo($motif_length,2,2)_motif/XML_of_DeNovo($motif_length,2,2)_motif.xml"
    id="$tmp_dir/$tag.fa"
    f=1.0
    outdir=$tmp_dir
    bs=false`);
    for line in eachline(string(tmp_dir,"/Motif_hits_from_SequenceScan(1.0).BED"))
        container = push!(container, parse(Float64, split(line)[5]))
    end
    run(`rm $tmp_dir/Binding_sites_from_SequenceScan\(1.0\).txt $tmp_dir/Motif_hits_from_SequenceScan\(1.0\).BED $tmp_dir/protocol_scan.txt $tmp_dir/$tag.fa`);
    return(log10.(container))
end


function make_inmode(path_to_inmode::String, path_to_java::String, motif_length::Int64, order::Int64, tmp_dir::String)
    run(`$path_to_java -Xmx4096m -Xms1024m --add-modules java.xml.bind -jar $path_to_inmode denovo
    i="$tmp_dir/train.fa"
    m=$motif_length
    outdir=$tmp_dir
    mo=$order`)
end


function bootstrap_inmode(sites::Array{Array{Char, 1}, 1}, path_to_inmode::String, path_to_java::String, tmp_dir::String, order::Int64, size_of::Int64)
    true_scores = Float64[]
    false_scores = Float64[]
    number_of_sites = length(sites)
    motif_length = length(sites[1])

    for i in 1:10

        if !isdir(tmp_dir)
            mkdir(tmp_dir)
        end

        index_train = Random.randsubseq(1:number_of_sites, 0.9)
        index_test = setdiff(1:number_of_sites, index_train)
        index_shuffle = Random.rand(1:number_of_sites, size_of)

        write_fasta(sites[index_train], tmp_dir, "train")
        write_fasta(sites[index_test], tmp_dir, "test")
        write_fasta(Random.shuffle.(sites[index_shuffle]), tmp_dir, "shuffled")

        make_inmode(path_to_inmode, path_to_java, motif_length, order, tmp_dir)
        true_scores = vcat(true_scores, calculate_scores(path_to_inmode, path_to_java, motif_length, tmp_dir, "test"))
        false_scores = vcat(false_scores, calculate_scores(path_to_inmode, path_to_java, motif_length, tmp_dir, "shuffled"))

        #rm(tmp_dir, recursive=true)

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
            help = "path to fasta file with sites"
            required = true
        "inmode"
            help = "path to inmode"
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
        "--java", "-j"
            help = "path to java"
            required = false
            arg_type = String
            default = "java"
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
    path_to_inmode = args["inmode"]
    path_to_java = args["java"]

    sites = read_sites(path)
    df = bootstrap_inmode(sites, path_to_inmode, path_to_java, tmp_dir, order, size_of)
    CSV.write(out, df, delim='\t')

end

main()
