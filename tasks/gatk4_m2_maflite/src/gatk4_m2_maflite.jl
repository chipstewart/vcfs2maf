#!/usr/local/bin/julia
# ARGS=["REBC-AC8L-TP","REBC-AC8L-NT","tmp1.tsv","REBC-AC8L-TP-NT.M2_maflite.tsv"]
# ARGS=["SU2CLC-MGH-1047-TM-01","SU2CLC-MGH-1047-BL-01","tmp3.tsv","RP-1066_SU2CLC-MGH-1047"]
# ARGS=["REBC-AF8C-NT1-A-1-1-D-A649-36","REBC-AF8C-TTP1-A-1-1-D-A649-36","tmp3.tsv","REBC-AF8C-NT-TP","37"]

using DataFrames, CSV, Statistics, DelimitedFiles
tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
pair_id=ARGS[4]
build_id=ARGS[5]

file1a=pair_id*".raw.tsv"

isString(x::Number)=false
#isString(x::Array{>:Missing}.NAtype)=false
isString(x::AbstractString)=true

#df = readtable(file1)
df = CSV.read(file1; delim='\t')





# size(df)
# describe(df)
# print(df)
FIELDS=names(df)
FIELDX=[:x,:QUAL,:ID,:NORMAL_FT,:NORMAL_MMQ,:NORMAL_PID, :NORMAL_MBQ, :NORMAL_GQ, :NORMAL_OBF, :NORMAL_OBQRC,:NORMAL_AF,:NORMAL_FOXOG,
 :NORMAL_SA_MAP_AF, :NORMAL_MPOS, :NORMAL_PL, :NORMAL_OBAM, :NORMAL_OBQ, :NORMAL_ALT_F2R1, :NORMAL_SA_POST_PROB, :NORMAL_REF_F2R1,
 :NORMAL_GT, :NORMAL_MCL, :NORMAL_MFRL, :NORMAL_OBP, :NORMAL_OBAMRC, :NORMAL_PGT, :NORMAL_REF_F1R2,  :NORMAL_ALT_F1R2,
 :TUMOR_FT, :TUMOR_PID,  :TUMOR_GQ, :TUMOR_OBF, :TUMOR_OBQRC, :TUMOR_SA_MAP_AF, :TUMOR_PL, :TUMOR_OBAM, :NORMAL_SB, 
 :TUMOR_OBQ, :TUMOR_SA_POST_PROB,:TUMOR_GT,:TUMOR_OBP, :TUMOR_OBAMRCM, :TUMOR_PGT, :TUMOR_SB, :OCM, :MMQ, :MBQ, :SEQQ,:STRQ, :MPOS, :RPA,:GERMQ, 
 :MFRL, :CONTQ, :ECNT,:NALOD, :STRANDQ, :TUMOR_PS, :NORMAL_PS, :NORMAL_F2R1_REF, :NORMAL_F2R1_ALT, :NORMAL_F1R2_REF, :NORMAL_F1R2_ALT]

for i=1:length(FIELDX)
    if (FIELDX[i] in FIELDS)
        deletecols!(df, FIELDX[i])
    end
end



for c in names(df)
    println(c)
    k=map(x->ismissing(x),df[c])
    if mean(k)==1.0
        deletecols!(df, c)
        continue
    end
    if isa(df[c],Array{Union{Missing, String},1}) 
        df[k,c] = ""
    end
end
#rename!(df, [:normal_tumor_alt_count, :normal_tumor_ref_count], [:t_alt_count, :t_ref_count])
head(df)
a= df[:CHRO]
if !isa(a[1],Int)
    a=replace(a,"X"=>"23","Y"=>"24","MT"=>"25","M"=>"25")
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
# print(df)
sort!(df, [:a, :POS])
deletecols!(df, [:a])

println("")
println(file1a)
println("")

#for c in [:TUMOR_ALT_F1R2,:TUMOR_ALT_F2R1,:TUMOR_REF_F1R2,:TUMOR_REF_F2R1,:TUMOR_FOXOG]
for c in [:TUMOR_F1R2_ALT,:TUMOR_F2R1_ALT,:TUMOR_F1R2_REF,:TUMOR_F2R1_REF]
    if ~isString(df[1,c])
        df[c] = map(x -> string(x),df[c])
    end
    k=findall(map(x->x=="NA",df[c]))
    df[k,c] = ""
end


open(file1a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Matrix,df), '\t')
end

a=df[:ALT]
r=df[:REF]
n=length(a)
ar=fill("",n)
INDEL=fill(0,n)
SNV=fill(0,n)
for i in 1:n
    #println(i)
    ar[i]=string(a[i],":",r[i])
    INDEL[i]=length(a[i])!=length(r[i])
    SNV[i]=length(a[i])==length(r[i])
end

#l1=map(x -> length(x), ar)
#df[:INDEL]=map(x -> x>3,l1)
df[:INDEL]=INDEL.>0
df[:Variant_Type]=fill("SNV",n)
df[df[:INDEL],:Variant_Type]="INDEL"
dfindel=df[df[:INDEL],:]
dfsnv=df[SNV.>0,:]
#deletecols!(dfsnv, [:INDEL])
#delete!(dfindel, [:INDEL])
#delete!(df, [:INDEL])
df2=df

# set ALT0 to original ALT,shorten ALT to at most 5 bases, then copy back after merge

df2[:ALT0]=df2[:ALT]
alt0=df2[:ALT0]
alt=fill("",n)
n=length(alt0)
l0=map(x -> length(x), alt0)
for i in 1:n
    l1=min(l0[i],5)
    alt1=string(alt0[i])
    alt[i]=alt1[1:l1]
end
df2[:ALT]=alt



#delete!(df, [:INDEL, :Variant_Type,:ALT0])



# label M2 with M2 flag
a= df[:CHRO]
if !isa(a[1],Int)
    #a=map(x -> replace(x,r"[X]", "23"), a)
    #a=map(x -> replace(x,r"[Y]", "24"), a)
    #a=map(x -> replace(x,r"[MT]", "25"), a)
    #a=map(x -> replace(x,r"[M]", "25"), a)
    a=replace(a,"X"=>"23","Y"=>"24","MT"=>"25","M"=>"25")  
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
#print(df)
sort!(df, [:a, :POS])
deletecols!(df, [:a])



for c in names(df)
     println(c)
    #k=ismissing(df[c])
    k=map(x->ismissing(x),df[c])
    if mean(k)==1.0
        delete!(df, c)
        continue
    end
    println(typeof(df[1,c]))
    if isString(df[1,c])
        df[k,c] = ""
    end
end

# multi-alt allele events 
k=findall(map(x-> length(split(x,':'))>1, df[:ALT]))
df[k,:ALT]=map(x-> split(x,':')[1], df[k,:ALT])
df[k,:POPAF]=map(x-> split(x,':')[1], df[k,:POPAF])
df[k,:NLOD]=map(x-> split(x,':')[1], df[k,:NLOD])
df[k,:TLOD]=map(x-> split(x,':')[1], df[k,:TLOD])
df[k,:TUMOR_AF]=map(x-> split(x,':')[1], df[k,:TUMOR_AF])



df[:end]=df[:POS]


# start for indels is POS+1
kdel = map(x-> length(x)>1, df[:REF]).*df[:INDEL]
p=df[:POS]
p=p+kdel
df[:POS]=p
k = findall(kdel)
ref=df[k,:REF]
ref=map(x->x[2:end],ref)
alt=df[k,:ALT]
alt=map(x->x[2:end],alt)
k1 = findall(map(x-> length(x)<1, alt))
alt[k1]=map(x-> "-", alt[k1])
df[k,:REF]=ref
df[k,:ALT]=alt
dref = map(x-> length(x), ref)
dalt = map(x-> length(x), alt)
ddel = dref
df[k,:end]=df[k,:POS]+ddel




kins = map(x-> length(x)>1, df[:ALT]).*df[:INDEL]
k = findall(kins)
ref=df[k,:REF]
ref=map(x->x[2:end],ref)
k1 = findall(map(x-> length(x)<1, ref))
ref[k1]=map(x-> "-", ref[k1])
alt=df[k,:ALT]
alt=map(x->x[2:end],alt)
df[k,:REF]=ref
df[k,:ALT]=alt
df[k,:end]=df[k,:POS].+1



kmnp = map(x-> length(x)>1, df[:ALT]).*map(x->!x,df[:INDEL])
k = findall(kmnp)
ref=df[k,:REF]
dref = map(x-> length(x), ref)
df[k,:end]=df[k,:POS]+dref.-1



# maflite fields 
# SNP: build    chr start   end ref_allele  tum_allele1 tum_allele2 tumor_barcode   normal_barcode  tumor_f init_t_lod  t_lod_fstar t_alt_count t_ref_count judgement
# INDEL: build  contig  start_position  end_position    ref_allele  alt_allele1 alt_allele2 tumor_name  normal_name tumor_f n_ref_count n_alt_count t_ref_count t_alt_count init_n_lod  init_t_lod  judgement


# add build column
df[:build]=fill(build_id,size(df,1))
df[:tumor_barcode]=fill(tumor_id,size(df,1))
df[:normal_barcode]=fill(normal_id,size(df,1))
df[:judgement]=fill("KEEP",size(df,1))


# allele counts
df[:n_alt_count]=df[:NORMAL_AD_ALT]
df[:n_ref_count]=df[:NORMAL_AD_REF]
#$df[:t_lod_fstar]=df[:TLOD]
df[:tumor_f]=df[:TUMOR_AF]
df[:t_alt_count]=df[:TUMOR_AD_ALT]
df[:t_ref_count]=df[:TUMOR_AD_REF]
rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])#from = [:DB, :HCNT, :MAX_ED,:RPA,:STR, :RU, :NLOD,:TUMOR_ALT_F1R2,:MIN_ED,:ECNT, :TLOD, :TUMOR_PID,:TUMOR_AF,:TUMOR_FOXOG, :TUMOR_QSS, :TUMOR_ALT_F2R1,:TUMOR_AD_REF,:TUMOR_AD_ALT, :TUMOR_REF_F2R1, :TUMOR_REF_F1R2,:NORMAL_AF,:NORMAL_AD_REF,:NORMAL_AD_ALT]
#from = ["DB","HCNT","MAX_ED","RPA","STR","RU","NLOD","TUMOR_ALT_F1R2","MIN_ED","ECNT","TLOD","TUMOR_PID","TUMOR_AF","TUMOR_FOXOG","TUMOR_QSS","TUMOR_ALT_F2R1","TUMOR_AD_REF","TUMOR_AD_ALT","TUMOR_REF_F2R1","TUMOR_REF_F1R2","NORMAL_AF","NORMAL_AD_REF","NORMAL_AD_ALT"]
#println(from)

#for c in names(df)
#    c1 = string(c)
#    println(c1)
#    if length(find(Bool[contains(c1,i) for i in from]))>0
#        m2c = Symbol(string("M2_",c1))
#        rename!(df,c,m2c)
#    end
#end
FIELDS=names(df)
FIELDX=[:NORMAL_DP, :NORMAL_AD_REF, :NORMAL_AD_ALT, :INDEL, :Variant_Type, :ALT0, 
 :TUMOR_DP, :TUMOR_AF, :TUMOR_AD_REF, :TUMOR_AD_ALT, :DP]
for i=1:length(FIELDX)
    if (FIELDX[i] in FIELDS)
        deletecols!(df, FIELDX[i])
    end
end

maf=df
for c in names(maf)
    #println(c)
    if ~isString(maf[1,c])
        maf[c] = map(x -> string(x),maf[c])
    end
end


for c in names(maf)
    println(c)
    k=map(x->ismissing(x),df[c]) 
    if mean(k)==1.0
        delete!(maf, c)
        continue
    end
    k=findall(map(x->x=="NA",maf[c]))
    maf[k,c] = ""
    k=findall(map(x->uppercase(x)=="TRUE",maf[c]))
    maf[k,c] = "1"
    k=findall(map(x->uppercase(x)=="FALSE",maf[c]))
    maf[k,c] = "0"
end

maflite_all=pair_id*".m2.all.maflite.tsv"
open(maflite_all, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Matrix,maf), '\t')
end

kpass = find(map(x-> uppercase(x)=="PASS", maf[:FILTER]))
maf=maf[kpass,:]
delete!(maf,[:FILTER])


maflite_pass=pair_id*".m2.pass.maflite.tsv"
open(maflite_pass, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Matrix,maf), '\t')
end


