using CodecZlib
using DataFrames
using BetaBinomialModels
using Distributions
using Statistics
using Printf
using CSV

pa = "/nfs/turbo/lsa-bis/data_for_Kerby/illumina"
dp = "/nfs/turbo/lsa-bis/data_for_Kerby/illumina_fuse"
fn = "Sample_XXXX_CpG.bedGraph"

include("utils.jl")

invlogit(x) = 1 / (1 + exp(-x))

function fit_model(meth, total)
    bb = fit(BetaBinomSeq, meth, total; maxiter=10, fuse_mean_wt=1.0, fuse_icc_wt=1.0,
             step0=0.01, rho0=0.5, rhofac=1.1, rhomax=2.0, verbose=true)
    return bb
end

function fitstats(bb)
    rr = diff(findall(diff(bb.logit_mean) .!= 0))
    qmean = median(rr)
    rmean = mean(rr)
    rr = diff(findall(diff(bb.logit_icc) .!= 0))
    qicc = median(rr)
    ricc = mean(rr)
    return (median_mean_gap=qmean, mean_mean_gap=rmean, median_icc_gap=qicc, mean_icc_gap=ricc)
end

function main()
    chrom = "1"
    out = open(dp, "$(chrom).log", "w")
    for id in ["5554", "5555", "5663"]
        da = prep_data(id, chrom)
        bb = fit_model(da[:, :meth], da[:, :total])
        st = fitstats(bb)
        write(out, "$(id)\n")
        write(out, string(st))
        write(out, "\n")
        n = length(da[:, :meth])
        logit_mean, logit_icc = coef(bb)
        rslt = DataFrame(chrom=fill(chrom, n), startpos=da[:, :startpos], endpos=da[:, :endpos],
                         mean=invlogit.(logit_mean), icc=invlogit.(logit_icc),
                         meth=da[:, :meth], total=da[:, :total])
        open(GzipCompressorStream, joinpath(dp, "$(id)_$(chrom).csv.gz"), "w") do io
            CSV.write(io, rslt)
        end
    end
    close(out)
end

main()
