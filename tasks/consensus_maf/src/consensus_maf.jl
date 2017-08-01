#!/usr/local/bin/julia
# ARGS=["REBC-AC8L-TP","REBC-AC8L-NT","sample.mutect.maflite.txt","REBC-AC8L-TP-NT.m2_maflite.tsv","REBC-AC8L-TP-NT.Strelka_maflite.tsv","REBC-AC8L-TP-NT.SvABA_maflite.tsv","M1","M2","STRELKA","SVABA","REBC-AC8L-TP-NT.merged_maflite.tsv"]
using DataFrames
in_maf_file=ARGS[1]
NALG=ARGS[2]
LONG=ARGS[3]
out_maf_file=ARGS[4]

print(ARGS)

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

maflite_fields=["build","chr","start","end","ref_allele","tum_allele1","tum_allele2","tumor_barcode","normal_barcode","t_alt_count","t_ref_count","judgement","n_alt_count","n_ref_count","tumor_f"]
maflite_symbols=map(x->Symbol(x),maflite_fields)


if !isfile(file1)
    println("missing input maf ", in_maf_file, " .")
end
df1 = readtable(in_maf_file,separator ='\t')
if Symbol("_end") in names(df1)
	rename!(df1, [:_end], [:end])
end

event_length =  map(x -> length(x),df1[:Reference_Allele])+map(x -> length(x),df1[:Tumor_Seq_Allele2])-1
pass =  map(x -> x>=NALG,df1[:NALG]) | map(x -> x>=LONG),event_length)

for c in names(df1)
	if ~isString(df1[1,c])
		df1[c] = map(x -> string(x),df1[c])
	end
end

k=find(pass)
maf=df[k,:]

open(consensus_maf_file, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end

