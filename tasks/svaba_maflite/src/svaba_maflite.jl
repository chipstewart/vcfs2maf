 #!/usr/local/bin/julia
#ARGS
# ARGS=["THCA-EM-A2CN-TP","THCA-EM-A2CN-NB","THCA-EM-A2CN-TP-NB.SvABA.INDEL.tsv","THCA-EM-A2CN-TP-NB.maflite.tsv"]
# ARGS=["THCA-EM-THCA-BJ-A191-TP","THCA-BJ-A191-NB","THCA-BJ-A191-TP-NB.SvABA.INDEL.tsv","THCA-BJ-A191-TP-NB.maflite.tsv"]
# ARGS=["THCA-DJ-A13W-TP","THCA-DJ-A13W-NB","THCA-DJ-A13W-TP-NB.SvABA.INDEL.tsv","THCA-DJ-A13W-TP-NB.SvABA_maflite.tsv"]
# ARGS=["RP-1066_SU2CLC-MGH-1048-TM-01","RP-1066.SU2CLC-MGH-1048-BL-01","xxx/RP-1066.SU2CLC-MGH-1048-TM-01.bam", "xxx/RP-1066.SU2CLC-MGH-1048-BL-01.bam", "SU2CLC-MGH-1048_1.SvABA.INDEL.tsv","SU2CLC-MGH-1048_1.SvABA_maflite.tsv","2"]

using DataFrames

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

tumor_id=ARGS[1]
normal_id=ARGS[2]
tumor_bam=ARGS[3]
tumor_bam=basename(tumor_bam)
normal_bam=ARGS[4]
normal_bam=basename(normal_bam)
file1=ARGS[5]
maflite=ARGS[6]

max_normal_alt_count=ARGS[7]
if isString(max_normal_alt_count)
    max_normal_alt_count=parse(Int64,max_normal_alt_count)
end

df = readtable(file1,normalizenames=false)
# size(df)
# describe(df)
delete!(df, [:FILTER,:QUAL,:ID,:SCTG])
#sampleN_NALT,:sampleN_NALT_SR,:sampleN_NREF,:sampleN_NALT_RP,:sampleN_READ_ID,:sampleT_NALT_SR,:sampleT_NREF,:sampleT_NALT_RP,:sampleT_NALT,:sampleT_READ_ID,:x])
# print(df)


for c in names(df)
    println(c)
    c1=string(c)
    k1=searchindex(c1,tumor_bam)
    if k1>0
        kk=search(c1,tumor_bam)
        k2=1+kk[endof(kk)]
        k3=length(c1)
        c2=string("tumor_",c1[k2+1:k3])
        c3=Symbol(c2)
        rename!(df,c, c3)
    end
    k1=searchindex(c1,normal_bam)
    if k1>0
        kk=search(c1,normal_bam)
        k2=1+kk[endof(kk)]
        k3=length(c1)
        c2=string("normal_",c1[k2+1:k3])
        c3=Symbol(c2)
        rename!(df,c, c3)
    end    
end  

df=df[df[:normal_AD].<=max_normal_alt_count,:]

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
    a=map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),a)
    #a=map(x -> tryparse(Float64,x), a)
    #xchrom=map(x -> isnull(x), a)
    #a[find(xchrom)]=0.0.*find(xchrom)
    #a=DataArray(map(x -> get(x, NaN), a), map(isnull, a))
end
df[:a] = a
df=df[df[:a].>0.1,:]

# print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])

# SvABA indel POS is first base of ALT (*NOT* VCF standard) 
a0=df[:ALT]
r0=df[:REF]
p=df[:POS]
a = map(x-> x[2:end], a0)
r = map(x-> x[2:end], r0)
lr = map(x-> length(x), r)
la = map(x-> length(x), a)
# insertions
kr0=find(map(x-> x==0,lr))
r[kr0]=fill("-",length(kr0))
# deletions
ka0=find(map(x-> x==0,la))
a[ka0]=fill("-",length(ka0))
# insertions default
df[:end]=p+1
# deletions
[df[:POS] df[:end] df[:REF] df[:ALT] p la lr]
df[ka0,:POS]=p[ka0]+1
[df[:POS] df[:end] df[:REF] df[:ALT] p la lr]
# mysterious incremented p - compensate by subtracting 1 
df[ka0,:end]=p[ka0]+lr[ka0]-1 
[df[:POS] df[:end] df[:REF] df[:ALT] p la lr]

df[:REF]=r;
df[:ALT]=a;

# ref counts
nref=df[:normal_DP]-df[:normal_AD]
tref=df[:tumor_DP]-df[:tumor_AD]

df[:n_ref_count]=nref
df[:t_ref_count]=tref
rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])
rename!(df, [:tumor_AD,:normal_AD,], [:t_alt_count, :n_alt_count])
#rename!(df, [:REPSEQ,:NM, :MAPQ, :DBSNP,:LOD,:SPAN],  [:SvABA_REPSEQ, :SvABA_NM,:SvABA_MAPQ,:SvABA_DBSNP,:SvABA_LOD,:SvABA_SPAN])


for c in names(df)
    println(c)
    c1=string(c)
    if searchindex(c1,"normal_PL")>0
        delete!(df, [c])
    end
    if searchindex(c1,"normal_GQ")>0
        delete!(df, [c])
    end
    if searchindex(c1,"normal_GT")>0
        delete!(df, [c])
    end
    if searchindex(c1,"normal_LR")>0
        delete!(df, [c])
    end
    if searchindex(c1,"tumor_PL")>0
        delete!(df, [c])
    end
    if searchindex(c1,"tumor_GQ")>0
        delete!(df, [c])
    end
    if searchindex(c1,"tumor_GT")>0
        delete!(df, [c])
    end
    if searchindex(c1,"tumor_LR")>0
        delete!(df, [c])
    end
end

#for c in names(df)
#    println(c)
#    c1=string(c)
#    if searchindex(c1,"tumor_") == 1
#        c2=string("SvABA_",c1)
#        rename!(df,c,Symbol(c2))
#    end
#    if searchindex(c1,"normal_") == 1
#        c2=string("SvABA_",c1)
#        rename!(df,c,Symbol(c2))
#    end
#end

[df[:start] df[:end] df[:ref_allele] df[:alt_allele]]



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
