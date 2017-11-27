#!/usr/local/bin/julia
# ARGS=["REBC-AC97-TP","REBC-AC97-NB","REBC-AC97-TP-NB.mutect.maflite.txt","REBC-AC97-TP-NT.m1_maflite.tsv"]

using DataFrames
tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
maflite=ARGS[4]

df = readtable(file1,separator='\t')
# size(df)
# describe(df)
# print(df)
delete!(df, [:tum_allele1])
rename!(df, [:_end, :tum_allele2], [:end, :alt_allele])

for c in names(df)
    println(c)
    k=isna(df[c])
    if mean(k)==1.0
		delete!(df, c)
		continue
	end
    if isa(df[c],DataArray{String,1}) 
        df[k,c] = ""
    end
end


isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true



maf=df
for c in names(maf)
    #println(c)
    if ~isString(maf[1,c])
        maf[c] = map(x -> string(x),maf[c])
    end
end

maf[:tumor_barcode]=tumor_id
maf[:normal_barcode]=normal_id

for c in names(maf)
    println(c)
    k=isna(maf[c])
    if mean(k)==1.0
        delete!(maf, c)
        continue
    end
    k=find(map(x->x=="NA",maf[c]))
    maf[k,c] = ""
    k=find(map(x->uppercase(x)=="TRUE",maf[c]))
    maf[k,c] = "1"
    k=find(map(x->uppercase(x)=="FALSE",maf[c]))
    maf[k,c] = "0"
end


open(maflite, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
