 #!/usr/local/bin/julia
#ARGS
using DataFrames
# tumor_id="T"
# normal_id="N"
# file1="/Volumes/64GB_2015/Work/REBC/test1.tsv"
# file1a="/Volumes/64GB_2015/Work/REBC/test1a.tsv"
# file2="/Volumes/64GB_2015/Work/REBC/test2.tsv"
# file2a="/Volumes/64GB_2015/Work/REBC/test2a.tsv"
# maflite="/Volumes/64GB_2015/Work/REBC/test_M1_M2.tsv"
tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
file2=ARGS[4]
maflite=ARGS[5]
file1a="vcf_1.tsv"
file2a="vcf_2.tsv"

df = readtable(file1)
#size(df)
#describe(df)
#print(df)
#head(df)

rename!(df, [:normal_tumor_alt_count, :normal_tumor_ref_count], [:t_alt_count, :t_ref_count])
delete!(df, [:FILTER, :judgement,:x,:normal_tumor_GT,:QUAL])
#head(df)
#a=df[:CHRO]
a = df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
#print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])
#writetable("/Volumes/64GB_2015/Work/REBC/test1a.tsv", df1, header = true,quotemark='')
open(file1a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Array,df), '\t')
end

df1=df

df = readtable(file2)
# size(df)
# describe(df)
# print(df)
delete!(df, [:FILTER,:QUAL, :NORMAL_PID,:x,:TUMOR_GT,:TUMOR_PL,:TUMOR_GQ,:TUMOR_DP,:TUMOR_PGT,:NORMAL_FOXOG,:NORMAL_QSS,:NORMAL_ALT_F2R1,:NORMAL_ALT_F1R2,:NORMAL_REF_F2R1,:NORMAL_REF_F1R2,:NORMAL_GT,:NORMAL_PGT,:NORMAL_GQ,:NORMAL_PL,:NORMAL_DP])
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
#rename!(df, [:normal_tumor_alt_count, :normal_tumor_ref_count], [:t_alt_count, :t_ref_count])
head(df)
a= df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
# print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])
open(file2a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Array,df), '\t')
end

df2=df

df= join(df1, df2, on =  [:CHRO, :POS, :REF, :ALT], kind = :outer)

# label M1 with M1 flag
df[:M1]=~isna(df[:n_alt_count])
df[:M2]=~isna(df[:TLOD])
a= df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
#print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

for c in names(df)
     println(c)
    k=isna(df[c])
    if mean(k)==1.0
		delete!(df, c)
		continue
	end
    println(typeof(df[1,c]))
    if isString(df[1,c])
    	df[k,c] = ""
    end
end

# maflite fields 
# SNP: build	chr	start	end	ref_allele	tum_allele1	tum_allele2	tumor_barcode	normal_barcode	tumor_f	init_t_lod	t_lod_fstar	t_alt_count	t_ref_count	judgement
# INDEL: build	contig	start_position	end_position	ref_allele	alt_allele1	alt_allele2	tumor_name	normal_name	tumor_f	n_ref_count	n_alt_count	t_ref_count	t_alt_count	init_n_lod	init_t_lod	judgement

df[:end]=df[:POS]

# start for indels is POS+1
kdel = map(x-> length(x)>1, df[:REF])
p=df[:POS]
p=p+kdel
df[:POS]=p
k = find(kdel)
ref=df[k,:REF]
ref=map(x->x[2:end],ref)
alt=df[k,:ALT]
alt=map(x->x[2:end],alt)
k1 = find(map(x-> length(x)<1, alt))
alt[k1]=map(x-> "-", alt[k1])
df[k,:REF]=ref
df[k,:ALT]=alt
dref = map(x-> length(x), ref)
dalt = map(x-> length(x), alt)
ddel = dref
df[k,:end]=df[k,:POS]+ddel




kins = map(x-> length(x)>1, df[:ALT])
# p=df[:POS]
# p=p+kins
# df[:POS]=p
# df[find(kins),:POS]
k = find(kins)
ref=df[k,:REF]
ref=map(x->x[2:end],ref)
k1 = find(map(x-> length(x)<1, ref))
ref[k1]=map(x-> "-", ref[k1])
alt=df[k,:ALT]
alt=map(x->x[2:end],alt)
df[k,:REF]=ref
df[k,:ALT]=alt
df[k,:end]=df[k,:POS]+1



# add build column
df[:build]=fill("37",size(df,1))
df[:tumor_barcode]=fill(tumor_id,size(df,1))
df[:normal_barcode]=fill(normal_id,size(df,1))
df[:judgement]=fill("KEEP",size(df,1))

# M1 counts
k2=find(isna(df[:t_alt_count]))
df[k2,:ID]=df[k2,:ID_1]
df[k2,:n_alt_count]=df[k2,:NORMAL_AD_ALT]
df[k2,:n_ref_count]=df[k2,:NORMAL_AD_REF]
df[k2,:t_lod_fstar]=df[k2,:TLOD]
df[k2,:tumor_f]=df[k2,:TUMOR_AF]
df[k2,:init_t_lod]=df[k2,:NORMAL_AD_ALT]
df[k2,:t_alt_count]=df[k2,:TUMOR_AD_ALT]
df[k2,:t_ref_count]=df[k2,:TUMOR_AD_REF]
rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])
#from = [:ID_1,:DB, :HCNT, :MAX_ED,:RPA,:STR, :RU, :NLOD,:TUMOR_ALT_F1R2,:MIN_ED,:ECNT, :TLOD, :TUMOR_PID,:TUMOR_AF,:TUMOR_FOXOG, :TUMOR_QSS, :TUMOR_ALT_F2R1,:TUMOR_AD_REF,:TUMOR_AD_ALT, :TUMOR_REF_F2R1, :TUMOR_REF_F1R2,:NORMAL_AF,:NORMAL_AD_REF,:NORMAL_AD_ALT]
from = ["ID_1","DB","HCNT","MAX_ED","RPA","STR","RU","NLOD","TUMOR_ALT_F1R2","MIN_ED","ECNT","TLOD","TUMOR_PID","TUMOR_AF","TUMOR_FOXOG","TUMOR_QSS","TUMOR_ALT_F2R1","TUMOR_AD_REF","TUMOR_AD_ALT","TUMOR_REF_F2R1","TUMOR_REF_F1R2","NORMAL_AF","NORMAL_AD_REF","NORMAL_AD_ALT"]
println(from)

for c in names(df)
    c1 = string(c)
    println(c1)
    if length(find(Bool[contains(c1,i) for i in from]))>0
        m2c = Symbol(string("M2_",c1))
        rename!(df,c,m2c)
    end
end

#rename!(df, [:ID_1,:DB, :HCNT, :MAX_ED,:RPA,:STR, :RU, :NLOD,:TUMOR_ALT_F1R2], [:M2_ID,:M2_DB, :M2_HCNT, :M2_MAX_ED,:M2_RPA,:M2_STR, :M2_RU, :M2_NLOD,:M2_TUMOR_ALT_F1R2])
#rename!(df, [:MIN_ED,:ECNT, :TLOD, :TUMOR_PID,:TUMOR_AF,:TUMOR_FOXOG, :TUMOR_QSS, :TUMOR_ALT_F2R1], [:M2_MIN_ED,:M2_ECNT, :M2_TLOD, :M2_TUMOR_PID,:M2_TUMOR_AF,:M2_TUMOR_FOXOG, :M2_TUMOR_QSS, :M2_TUMOR_ALT_F2R1])
#rename!(df, [:TUMOR_AD_REF,:TUMOR_AD_ALT, :TUMOR_REF_F2R1, :TUMOR_REF_F1R2,:NORMAL_AF,:NORMAL_AD_REF,:NORMAL_AD_ALT],[:M2_TUMOR_AD_REF,:M2_TUMOR_AD_ALT, :M2_TUMOR_REF_F2R1, :M2_TUMOR_REF_F1R2,:M2_NORMAL_AF,:M2_NORMAL_AD_REF,:M2_NORMAL_AD_ALT])

#maf=df(:[])

#writetable("/Volumes/64GB_2015/Work/REBC/test_M1_M2.tsv", df, header = true,quotemark='"')

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
    k=find(map(x->uppercase(x)=="TRUE",maf[c]))
    maf[k,c] = "1"
    k=find(map(x->uppercase(x)=="FALSE",maf[c]))
    maf[k,c] = "0"
end


open(maflite, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
