#!/usr/local/bin/julia
# ARGS=["REBC-A90D-TP-NB-BI.intersection.maf","2","50","REBC-A90D-TP-NB-BI.consensus.maf"]
# ARGS=["THCA-BJ-A45K-TP-NB.filter.intersection.maf","2","50","THCA-BJ-A45K-TP-NB.consensus.maf"]
# ARGS=["/opt/test/REBC-A90K-TP-NB-BI.filter.intersection.maf","2","50","/opt/test/REBC-A90K-TP-NB-BI.consensus.maf"]
using DataFrames
in_maf_file=ARGS[1]
NALG=ARGS[2]
LONG=ARGS[3]
out_maf_file=ARGS[4]

print(ARGS)

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

if isString(NALG)
    NALG = parse(Float64,NALG)
end
if isString(LONG)
    LONG = parse(Float64,LONG)
end

if !isfile(in_maf_file)
    println("missing input maf ", in_maf_file, " .")
end

df1 = readtable(in_maf_file,separator ='\t')

nalgx=df1[:NALG]
nalg=fill(0,length(nalgx))

if all(map(x->isString(x),nalgx)) 
    ni=0
    for n1 in nalgx
        println(ni,n1)
        ni+=1
        nalg[ni]=maximum(map(x->parse(Int32,x),split(n1,"|";limit=5)))        
    end
else
    nalg=nalgx
end

delete!(df1, [:NALG])
df1[:NALG]=nalg

event_length =  map(x -> length(x),df1[:Reference_Allele])+map(x -> length(x),df1[:Tumor_Seq_Allele2])-1
pass =  map(x -> x>=NALG,nalg) | map(x -> x>=LONG,event_length)

for c in names(df1)
    if ~all(map(x->isString(x),df1[c]))
        df1[c] = map(x -> string(x),df1[c])
    end
end

for c in names(df1)
    k=find(map(x -> x=="NA",df1[c]))
    df1[k,c] = ""
end

k=find(pass)
maf=df1[k,:]

open(out_maf_file, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
