 #!/usr/local/bin/julia
#ARGS
# tumor_id="REBC-YQ-A90L-TP"; normal_id="REBC-A90L-NB"; file1="REBC-REBC-YQ-A90L.m1m2_maflite.tsv"; file2="REBC-REBC-YQ-A90L.SM.tsv";  maflite="/opt/test/REBC-A90L-TP-NT.m1m2_snowman_maflite.tsv"

using DataFrames
# tumor_id="REBC-A90I-TP"
# normal_id="REBC-A90I-NT"
# file1="/opt/test/REBC-A90I-TP-NT.m1m2_maflite.tsv"
# file2="/opt/test/REBC-A90I-TP-NT.broad-snowman.DATECODE.somatic.indel.tsv"
# maflite="/opt/test/REBC-A90I-TP-NT.m1m2_snowman_maflite.tsv"
# tumor_id="REBC-YQ-A90E"; normal_id="REBC-A90E-NT"; file1="REBC-REBC-YQ-A90E.m1m2_maflite.tsv"; file2="SM.tsv";  maflite="/opt/test/REBC-A90E-TP-NT.m1m2_snowman_maflite.tsv"
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

a = df[:chr]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    a=map(x -> replace(x,r"[MT]", "25"), a)
    a=map(x -> replace(x,r"[M]", "25"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
#print(df)
sort!(df, cols = [:a, :start])
delete!(df, [:a])
#writetable("/Volumes/64GB_2015/Work/REBC/test1a.tsv", df1, header = true,quotemark='')
rename!(df, [:_end], [:end])
df1=df

df = readtable(file2)
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

n=length(a)
df[:SM]=fill(1,n)

# ar=fill("",n)
# for i in 1:n
#     a1=a[i]
#     a1=a[2:end]
#     #println(i)
#     ar[i]=string(a[i],":",r[i])
# end

# l1=map(x -> length(x), ar)
#df[:INDEL]=map(x -> x>3,l1)
#df[:Variant_Type]=fill("SNV",n)
#df[df[:INDEL],:Variant_Type]="INDEL"

#df[:tumor_barcode]=fill(tumor_id,n)
#df[:normal_barcode]=fill(normal_id,n)
#delete!(dfsnv, [:INDEL])
#delete!(dfindel, [:INDEL])
#delete!(df, [:INDEL])

df[:n_ref_count]=df[:NCOV]-df[:NSPLIT]
df[:t_ref_count]=df[:TCOV]-df[:TSPLIT]
delete!(df, [:TCOV,:NCOV,:NCIGAR,:NFRAC])
rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])
rename!(df, [:TFRAC,:NSPLIT, :TSPLIT, :TCIGAR,:MAPQ,:BLACKLIST,:SPAN], [:tumor_f, :n_alt_count,:t_alt_count,:SM_TCIGAR,:SM_MAPQ,:SM_BLACKLIST,:SM_SPAN])

# df[:end]=df[:start]
# # start for indels is POS+1
# kdel = map(x-> length(x)>1, df[:ref_allele])
# p=df[:start]
# p=p+kdel
# df[:start]=p
# k = find(kdel)
# ref=df[k,:ref_allele]
# ref=map(x->x[2:end],ref)
# alt=df[k,:alt_allele]
# alt=map(x->x[2:end],alt)
# k1 = find(map(x-> length(x)<1, alt))
# alt[k1]=map(x-> "-", alt[k1])
# df[k,:ref_allele]=ref
# df[k,:alt_allele]=alt
# dref = map(x-> length(x), ref)
# dalt = map(x-> length(x), alt)
# ddel = dref
# df[k,:end]=df[k,:start]+ddel




# kins = map(x-> length(x)>1, df[:ref_allele])
# k = find(kins)
# ref=df[k,:ref_allele]
# ref=map(x->x[2:end],ref)
# k1 = find(map(x-> length(x)<1, ref))
# ref[k1]=map(x-> "-", ref[k1])
# alt=df[k,:alt_allele]
# alt=map(x->x[2:end],alt)
# df[k,:ref_allele]=ref
# df[k,:alt_allele]=alt
# df[k,:end]=df[k,:start]+1

df2=df

# set ALT0 to original ALT,shorten ALT to at most 5 bases, then copy back after merge

p=map(x ->string(x),df2[:start])
n=length(p)
for i=1:n
    p[i]=string(df2[i,:chr],":",p[i],"@",df2[i,:alt_allele])
end
df2[:key1]=p

p=map(x ->string(x),df1[:start])
n=length(p)
for i=1:n
    p[i]=string(df1[i,:chr],":",p[i],"@",df1[i,:alt_allele])
end
df1[:key1]=p

# merge SNV                                        
df= join(df1, df2,  on =  [:key1,:chr], kind = :outer)


# df2[:ALT0]=df2[:alt_allele]
# alt0=df2[:ALT0]
# n=length(alt0)
# alt=fill("",n)
# l0=map(x -> length(x), alt0)
# for i in 1:n
#     l1=min(l0[i],5)
#     alt1=string(alt0[i])
#     alt[i]=alt1[1:l1]
# end
# df2[:alt_allele]=alt

# df1[:ALT0]=df1[:alt_allele]
# alt0=df1[:ALT0]
# n=length(alt0)
# alt=fill("",n)
# l0=map(x -> length(x), alt0)
# for i in 1:n
#     l1=min(l0[i],5)
#     alt1=string(alt0[i])
#     alt[i]=alt1[1:l1]
# end
# df1[:alt_allele]=alt

# # df1 indels 
# ki=find(map(x-> '-' in x, df1[:ref_allele])|map(x-> '-' in x, df1[:alt_allele]))
# dfi=df1[ki,:]
# # merge indels
# df= join(dfi, df2,  on =  [:chr, :start,:alt_allele], kind = :outer)

# # df1 snvs
# q=(map(x-> '-' in x, df1[:ref_allele])|map(x-> '-' in x, df1[:alt_allele]))
# ks=find((map(x-> !contains(x,"-"), df1[:ref_allele])) & map(x-> !contains(x,"-"), df1[:alt_allele]))
# dfs=df1[ks,:]

# # merge SNV
# df= join(dfi, df2,  on =  [:chr, :start,:alt_allele], kind = :outer)
# df= join(df, dfs,  on =  [:chr, :start,:alt_allele], kind = :outer)
# # indels
# #df= join(df, dfindel,  on =  [:CHRO, :POS ], kind = :outer)

# label M1 with M1 flag
k=find(isna(df[:M1]))
df[k,:M1]=fill(false,length(k))
k=find(isna(df[:M2]))
df[k,:M2]=fill(false,length(k))
k=find(isna(df[:SM]))
df[k,:SM]=fill(false,length(k))

#k=find(1-(df[:M2]| df[:M1]))

# algorithm
a=4*df[:SM]+2*df[:M2]+df[:M1]
k=find(map(x -> x==4,a))

df[k,:start]=df[k,:start_1]
df[k,:ref_allele]=df[k,:ref_allele_1]
df[k,:alt_allele]=df[k,:alt_allele_1]
#df[k,:ALT0]=df[k,:ALT0_1]
#df[k,:tumor_barcode]=df[k,:tumor_barcode_1]
#df[k,:normal_barcode]=df[k,:normal_barcode_1]
df[k,:n_ref_count]=df[k,:n_ref_count_1]
df[k,:t_ref_count]=df[k,:t_ref_count_1]
df[k,:t_alt_count]=df[k,:t_alt_count_1]
df[k,:n_alt_count]=df[k,:n_alt_count_1]
df[k,:tumor_f]=df[k,:tumor_f_1]
df[k,:end]=df[k,:end_1]
df[k,:ID]=fill(".",length(k))

# put ALT0 back
#df[:alt_allele]=df[:ALT0]
#delete!(df, [:INDEL, :Variant_Type,:ALT0,:ALT0_1,:tumor_barcode_1,:normal_barcode_1])
delete!(df, [:start_1,:ref_allele_1,:alt_allele_1,:key1])
delete!(df, [:n_ref_count_1,:t_ref_count_1,:tumor_f_1,:n_alt_count_1,:t_alt_count_1,:end_1])



# label M1 with M1 flag

a= df[:chr]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    a=map(x -> replace(x,r"[MT]", "25"), a)
    a=map(x -> replace(x,r"[M]", "25"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
#print(df)
sort!(df, cols = [:a, :start])
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






# add build column
df[:build]=fill("37",size(df,1))
df[:tumor_barcode]=fill(tumor_id,size(df,1))
df[:normal_barcode]=fill(normal_id,size(df,1))
df[:judgement]=fill("KEEP",size(df,1))


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
