#!/usr/local/bin/julia
#ARGS
# ARGS=["SU2CLC-MGH-1158_T1", "/opt/test/SU2CLC-MGH-1158_T1.consensus.maf",  "Protein_Change:HGVSp_Short;SwissProt_acc_Id:SWISSPROT", ".copy_fields.maf"]
# ARGS=["SU2CLC-MGH-1158_T1", "/opt/test/SU2CLC-MGH-1158_T1.consensus.maf",  "Protein_Change:HGVSp_Short;i_oxoG_Q:oxoG_Q", ".copy_fields.maf"]

using DataFrames

id=ARGS[1]
in_file=ARGS[2]
copy_fields=ARGS[3]
extension=ARGS[4]

# parse copy_fields
cfs=split(copy_fields, "*")
nf=length(cfs)
cf = Dict()
for cf1 in cfs
  fs=split(cf1, ":")
  cf[fs[1]] = fs[2]
end
 
f = open(in_file);
skip = -1
x="#"
while x[1]=='#'
     x = readline(f)
     skip += 1
end
close(f)

#maf
df = readtable(in_file,separator='\t',skipstart = skip, nastrings=[""])

#f = names(df)
f=map(x->string(x),names(df))
k=keys(cf)

k2=intersect(k,f)

if length(k2)<length(k)
   for k1 in k
     if ! (k1  in k2)
        print(k1," missing")
     end
   end
   error("missing fields")
end

for k1 in k
    f1=symbol(k1)
    f2=symbol(cf[k1])
    df[f2]=df[f1]
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

out_file=string(id,extension)
open(out_file, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
