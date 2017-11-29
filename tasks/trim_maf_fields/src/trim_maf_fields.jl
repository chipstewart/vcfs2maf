#!/usr/local/bin/julia
#ARGS
# ARGS=["THCA-EM-A2CP-TP-NB", "THCA-EM-A2CP-TP-NB.annotated.maf",  "fields_to_remove.txt", "THCA-EM-A2CP-TP-NB.annotated.trim.maf"]

using DataFrames

id=ARGS[1]
in_file=ARGS[2]
remove_file=ARGS[3]
out_file=ARGS[4]
#comment_char=ARGS[5]
#separator=ARGS[6]

f = open(remove_file);
remove = readlines(f)
close(f)

f = open(in_file);
skip = -1
x="#"
while x[1]=='#'
     x = readline(f)
     skip += 1
end
close(f)

#Strelka SNV vcf
df = readtable(in_file,separator='\t',skipstart = skip, nastrings=[""])
f = names(df)
r=map(x -> Symbol(chomp(x)),remove)

k=findin(r, f)
r=r[k]

delete!(df, r)



isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

maf=df
for c in names(maf)
    #println(c)
    if ~isString(maf[1,c])
        maf[c] = map(x -> string(x),maf[c])
    end
    k=find(map(x -> isna(x),maf[c]))
    maf[k,c]=""
end


for c in names(maf)
    println(c)

    k=find(map(x->x=="NA",maf[c]))
    maf[k,c] = ""
    k=find(map(x->uppercase(x)=="__UNKNOWN__",maf[c]))
    maf[k,c] = ""
end

open(out_file, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end


