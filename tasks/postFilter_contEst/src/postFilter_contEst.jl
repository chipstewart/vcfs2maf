#!/usr/local/bin/julia
# ARGS=["REBC-A90K-TP-NB-BI","/opt/test/REBC-A90K-TP-NB-BI.filter.intersection.maf","0.05","0.1"]
using DataFrames
using StatsFuns
using Distributions
ID=ARGS[1]
INPUT_MAF_FILE=ARGS[2]
CONTEST=ARGS[3]
PTHRESHOLD=ARGS[4]

print(ARGS)

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

if isString(CONTEST)
    CONTEST = parse(Float64,CONTEST)
end
if isString(PTHRESHOLD)
    PTHRESHOLD = parse(Float64,PTHRESHOLD)
end

if !isfile(INPUT_MAF_FILE)
    println("missing input maf ", INPUT_MAF_FILE, " .")
end

df1 = readtable(INPUT_MAF_FILE,separator ='\t')

t_alt_count=df1[:t_alt_count]
t_ref_count=df1[:t_ref_count]

n=length(t_alt_count)
p=fill(0.0,n)
pass=fill(true,n)
for i = 1:n
    alt = max(0,t_alt_count[i])
    ref = max(0,t_ref_count[i])
    d = Beta(alt+1, ref+1)
    p[i] = cdf(d, CONTEST)
    pass[i]=p[i]<PTHRESHOLD
end

df1[:contEst] = fill(CONTEST,n)
df1[:pContEst] = p
df1[:pContEst_PASS] = pass


for c in names(df1)
    if ~all(map(x->isString(x),df1[c]))
        df1[c] = map(x -> string(x),df1[c])
    end
end

maf=df1
all_maf_file = string(ID,".postFilter_contEst.all.maf")
open(all_maf_file, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end

k=find(pass)
maf=df1[k,:]

pass_maf_file = string(ID,".postFilter_contEst.pass.maf")
open(pass_maf_file, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
