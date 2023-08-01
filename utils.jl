using CodecZlib
using CSV
using DataFrames

# Path to methylation data
pa = "/nfs/turbo/lsa-bis/data_for_Kerby/illumina"

function read_data(fp, chrom)

    chr, numer, denom = String[], Int[], Int[]
    startpos, endpos = Int[], Int[]
    open(fp) do io
        readline(io) # header
        for line in eachline(io)
            x = split(line, "\t")
            c = x[1][4:end]
            if c != chrom
                continue
            end
            meth = tryparse(Int, x[5])
            nmeth = tryparse(Int, x[6])
            spos = tryparse(Int, x[2])
            epos = tryparse(Int, x[3])
            if !(isnothing(meth) || isnothing(nmeth) || isnothing(spos) || isnothing(epos))
                tot = meth + nmeth
                push!(chr, c)
                push!(numer, meth)
                push!(denom, tot)
                push!(startpos, spos)
                push!(endpos, epos)
            end
        end
    end

    return chr, numer, denom, startpos, endpos
end

# Read the results of fused BB model.
function read_fuse(id, chr)
    dp = "/nfs/turbo/lsa-bis/data_for_Kerby/illumina_fuse"
    fn = "XXXX_$(chr).csv.gz"
    startpos, endpos, mn, icc = [], [], [], []
    open(GzipDecompressorStream, joinpath(dp, replace(fn, "XXXX"=>id))) do io
        readline(io) # header
        for line in eachline(io)
            x = split(line, ",")
            if x[1] != "1"
                return startpos, endpos, mn, icc
            end
            push!(startpos, parse(Int, x[2]))
            push!(endpos, parse(Int, x[3]))
            push!(mn, parse(Float64, x[4]))
            push!(icc, parse(Float64, x[5]))
        end
    end
    return DataFrame(startpos=startpos, endpos=endpos, mean=mn, icc=icc)
end

function prep_data(id, chrom)

    fn = "Sample_XXXX_CpG.bedGraph"
    fp = replace(fn, "XXXX"=>id)
    chr, meth, total, startpos, endpos = read_data(joinpath(pa, fp), chrom)
    return DataFrame(startpos=startpos, endpos=endpos, meth=meth, total=total)
end
