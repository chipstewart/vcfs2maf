 #!/usr/local/bin/julia
#ARGS
# ARGS=["THCA-EM-A2CN-TP","THCA-EM-A2CN-NB","THCA-EM-A2CN-TP-NB.SvABA.INDEL.tsv","THCA-EM-A2CN-TP-NB.maflite.tsv"]
# ARGS=["THCA-EM-THCA-BJ-A191-TP","THCA-BJ-A191-NB","THCA-BJ-A191-TP-NB.SvABA.INDEL.tsv","THCA-BJ-A191-TP-NB.maflite.tsv"]
# ARGS=["THCA-DJ-A13W-TP","THCA-DJ-A13W-NB","THCA-DJ-A13W-TP-NB.SvABA.INDEL.tsv","THCA-DJ-A13W-TP-NB.SvABA_maflite.tsv"]

using DataFrames

file1a="vcf_1.tsv"

tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
maflite=ARGS[4]

df = readtable(file1)
# size(df)
# describe(df)
delete!(df, [:FILTER,:QUAL,:ID,:sampleN_NALT,:sampleN_NALT_SR,:sampleN_NREF,:sampleN_NALT_RP,:sampleN_READ_ID,:sampleT_NALT_SR,:sampleT_NREF,:sampleT_NALT_RP,:sampleT_NALT,:sampleT_READ_ID,:x])
# print(df)

for c in names(df)
    #println(c)
    k=isna(df[c])
    if mean(k)==1.0
		delete!(df, c)
		continue
	end
	df[k,c] = ""
end
head(df)
a= df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    a=map(x -> replace(x,r"[MT]", "25"), a)
    a=map(x -> replace(x,r"[M]", "25"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
# print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])

# Snowman indels all have preceding base in ALT and REF 
a0=df[:ALT]
r0=df[:REF]
p=df[:POS]
a = map(x-> x[2:end], a0)
r = map(x-> x[2:end], r0)
lr = map(x-> length(x), r)
la = map(x-> length(x), a)
kr0=find(map(x-> x==0,lr))
r[kr0]=fill("-",length(kr0))
ka0=find(map(x-> x==0,la))
a[ka0]=fill("-",length(ka0))
df[:end]=p+1
df[ka0,:POS]=p[ka0]+1
df[ka0,:end]=p[ka0]+1+lr[ka0]
df[:REF]=r;
df[:ALT]=a;


# we don't know the overlap of reads between TSPLIT and TCIGAR, so we just take the max 
nalt=max(df[:NSPLIT],df[:NCIGAR])
talt=max(df[:TSPLIT],df[:TCIGAR])
nref=max(df[:NCOV],nalt)-nalt
tref=max(df[:TCOV],talt)-talt

df[:n_alt_count]=nalt
df[:t_alt_count]=talt
df[:n_ref_count]=nref
df[:t_ref_count]=tref
rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])
rename!(df, [:TFRAC,:NSPLIT, :TSPLIT, :TCIGAR,:MAPQ,:BLACKLIST,:SPAN,:TCOV,:NCOV,:NCIGAR,:NFRAC],  [:SvABA_TFRAC, :SvABA_NSPLIT,:SvABA_TSPLIT,:SvABA_TCIGAR,:SvABA_MAPQ,:SvABA_BLACKLIST,:SvABA_SPAN,:SvABA_TCOV,:SvABA_NCOV,:SvABA_NCIGAR,:SvABA_NFRAC])

df[:SvABA_NFRAC]=df[:SvABA_NFRAC]/1.0
k=find(map(x->x<0,df[:SvABA_NFRAC]))
df[k,:SvABA_NFRAC]=NaN

df[:SvABA_TFRAC]=df[:SvABA_TFRAC]/1.0
k=find(map(x->x<0,df[:SvABA_TFRAC]))
df[k,:SvABA_TFRAC]=NaN

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


for c in names(maf)
    println(c)
    k=isna(maf[c])
    if mean(k)==1.0
        delete!(maf, c)
        continue
    end
    k=find(map(x->x=="NA",maf[c]))
    maf[k,c] = ""
    k=find(map(x->x=="NaN",maf[c]))
    maf[k,c] = ""
    k=find(map(x->uppercase(x)=="TRUE",maf[c]))
    maf[k,c] = "1"
    k=find(map(x->uppercase(x)=="FALSE",maf[c]))
    maf[k,c] = "0"
end

n=length(maf[:chr])

# add build column
maf[:build]=fill("37",size(maf,1))
maf[:tumor_barcode]=fill(tumor_id,size(maf,1))
maf[:normal_barcode]=fill(normal_id,size(maf,1))
maf[:judgement]=fill("KEEP",size(maf,1))

open(maflite, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
