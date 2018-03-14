#!/usr/local/bin/julia
# ARGS=["/opt/test/REBC-AF8C-TP-NT.postFilter_contEst.all.maf","2","M1:1;M2:1;Strelka1+Strelka2:1;Snowman+SvABA:1","REBC-AF8C-TP-NT"]
using DataFrames
in_maf_file=ARGS[1]
NALG=ARGS[2]
ALGS=ARGS[3]
id=ARGS[4]
out_maf_all=string(id,".consensus.all.maf")
out_maf_pass=string(id,".consensus.pass.maf")

print(ARGS)

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

type BLOCK1
    algs::Array{SubString{String}}
    nalg::Int64
    weight::Float64
end

if isString(NALG)
    NALG = parse(Float64,NALG)
end

B=split(ALGS ,r";") 
nb=length(B)
BLOCKS=BLOCK1[]
i=0
for b in B
  print(b,'\n')
  i+=1
  q=split(b ,r":")
  print(q[1],'\n')
  a0=String(q[1])
  a1=split(a0 ,r"\+")
  w1=parse(Float64,q[2])
  b1=BLOCK1(a1,length(a1),w1)
  push!(BLOCKS, b1 )
end


if !isfile(in_maf_file)
    println("missing input maf ", in_maf_file, " .")
end

df1 = readtable(in_maf_file,separator ='\t')
N=length(df1[:Chromosome])


vote=fill(0,N)
for b in BLOCKS
    algs=b.algs
    nalg=b.nalg
    weight=b.weight
    vote1=fill(false,N)
    for a in algs
        a1=Symbol(a)
        v2=map(x -> string(x),df1[a1])
        # fix oncotator formatting issues
        v3=map(x -> ismatch(r"1",x),v2)
        print(a),print(countmap(v3),'\n')
        vote1 = vote1 | v3
        df1[a1]=v3
    end
    votew=weight*vote1
    vote = vote + votew
end
df1[:vote]=vote


nalgx=df1[:NALG]
nalg=fill(0,length(nalgx))

if all(map(x->isString(x),nalgx)) 
    ni=0
    for n1 in nalgx
        println(ni,n1)
        ni+=1
        nalg[ni]=round(maximum(map(x->parse(Float64,x),split(n1,"|";limit=5))))     
    end
else
    nalg=nalgx
end

delete!(df1, [:NALG])
df1[:NALG]=nalg

event_length =  map(x -> length(x),df1[:Reference_Allele])+map(x -> length(x),df1[:Tumor_Seq_Allele2])-1
pass =  map(x -> x>=NALG,vote) 

for c in names(df1)
    if ~all(map(x->isString(x),df1[c]))
        df1[c] = map(x -> string(x),df1[c])
    end
end

for c in names(df1)
    k=find(map(x -> x=="NA",df1[c]))
    df1[k,c] = ""
end

open(out_maf_all, "w") do f
    writedlm(f, reshape(names(df1), 1, length(names(df1))), '\t')
    writedlm(f, convert(Array,df1), '\t')
end

k=find(pass)
maf=df1[k,:]

open(out_maf_pass, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
