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
file3=ARGS[5]
file4=ARGS[6]
maflite=ARGS[7]

tumor_id="T"
normal_id="N"
file1="tmp1.tsv"
file2="tmp2.tsv"
file3="tmp3.tsv"
file4="tmp4.tsv"
maflite="maflite.tsv"

file1a="vcf_1.tsv"
file2a="vcf_2.tsv"
file3a="vcf_3.tsv"
file4a="vcf_4.tsv"

# M1 vcf
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

# M2 vcf
df = readtable(file2)
# size(df)
# describe(df)
# print(df)
delete!(df, [:FILTER,:QUAL, :NORMAL_PID,:x,:TUMOR_GT,:TUMOR_PL,:TUMOR_GQ,:TUMOR_DP,:TUMOR_PGT,:NORMAL_FOXOG,:NORMAL_QSS,:NORMAL_ALT_F2R1,:NORMAL_ALT_F1R2,:NORMAL_REF_F2R1,:NORMAL_REF_F1R2,:NORMAL_GT,:NORMAL_PGT,:NORMAL_GQ,:NORMAL_PL,:NORMAL_DP])
# print(df)

if ( :TUMOR_FOXOG in names(df) )  # FOXOG messed up by VariantAnnotator
	delete!(df, [:TUMOR_FOXOG])
end

for c in names(df)
    #println(c)
    k=isna(df[c])
    if mean(k)==1.0
		delete!(df, c)
		continue
	end
	if isa(df[c],DataArray{String,1})
      df[k,c] = ""
	end	
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

#Strelka SNV vcf
df = readtable(file3)
# size(df)
# describe(df)
# print(df)
delete!(df, [:FILTER,:QUAL, :SOMATIC,:QSS_NT,:TQSS_NT,:SGT])
# print(df)

for c in names(df)
    #println(c)
    k=isna(df[c])
    if mean(k)==1.0
		delete!(df, c)
		continue
	end
	if isa(df[c],DataArray{String,1})
      df[k,c] = ""
	end	
end
#rename!(df, [:normal_tumor_alt_count, :normal_tumor_ref_count], [:t_alt_count, :t_ref_count])
head(df)
# keep only CHRO = [1-Y]
df=df[find(in.(df[:CHRO], [["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]])),:]

# sort by CHRO and position
a= df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    #a=map(x -> replace(x,r"[MT]", "25"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
# print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])


# Strelka 'tier' for QSS -
# hypothesis: TQSS must be strelka's favorite tier for each particular call 
nt=df[:TQSS]
n=length(nt)
AA=df[:ALT]
RA=df[:REF]

q=df[:NORMAL_AU]
N1= fill("0",n)
q = map(x -> split(x,":"), q)
for i = 1:n
	print(i,"\n")
	N1[i]=q[i][nt[i]]
end
NA=map(x -> parse(Int32,x), N1)

q=df[:NORMAL_CU]
N1 = fill("0",n)
q = map(x -> split(x,":"), q)
for i = 1:n
	print(i,"\n")
	N1[i]=q[i][nt[i]]
end
NC=map(x -> parse(Int32,x), N1)

q=df[:NORMAL_GU]
N1 = fill("0",n)
q = map(x -> split(x,":"), q)
for i = 1:n
	print(i,"\n")
	N1[i]=q[i][nt[i]]
end
NG=map(x -> parse(Int32,x), N1)

q=df[:NORMAL_TU]
N1 = fill("0",n)
q = map(x -> split(x,":"), q)
for i = 1:n
	print(i,"\n")
	N1[i]=q[i][nt[i]]
end
NT=map(x -> parse(Int32,x), N1)

q=df[:TUMOR_AU]
N1= fill("0",n)
q = map(x -> split(x,":"), q)
for i = 1:n
	print(i,"\n")
	N1[i]=q[i][nt[i]]
end
TA=map(x -> parse(Int32,x), N1)

q=df[:TUMOR_CU]
N1 = fill("0",n)
q = map(x -> split(x,":"), q)
for i = 1:n
	print(i,"\n")
	N1[i]=q[i][nt[i]]
end
TC=map(x -> parse(Int32,x), N1)

q=df[TUMOR_GU]
N1 = fill("0",n)
q = map(x -> split(x,":"), q)
for i = 1:n
	print(i,"\n")
	N1[i]=q[i][nt[i]]
end
TG=map(x -> parse(Int32,x), N1)

q=df[:TUMOR_TU]
N1 = fill("0",n)
q = map(x -> split(x,":"), q)
for i = 1:n
	print(i,"\n")
	N1[i]=q[i][nt[i]]
end
TT=map(x -> parse(Int32,x), N1)

df[:n_ref_count] = fill(0,n)
df[:n_alt_count] = fill(0,n)
df[:t_ref_count] = fill(0,n)
df[:t_alt_count] = fill(0,n)

NB=[NA NC NG NT]
TB=[TA TC TG TT]
b1 = ["A" "C" "G" "T"]
kr=find( in.(df[:REF],[b1]) )
mr=findfirst.([b1],df[:REF])
ka=find( in.(df[:ALT],[b1]) )
ma=findfirst.([b1],df[:ALT])

for i = 1:n
   df[i,:n_ref_count] = NB[i,mr[i]]
   df[i,:n_alt_count] = NB[i,ma[i]]
   df[i,:t_ref_count] = TB[i,mr[i]]
   df[i,:t_alt_count] = TB[i,ma[i]]
end
   
df[:QS] = df[:QSS]


open(file3a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Array,df), '\t')
end

df3=df

#Strelka INDEL vcf
df = readtable(file4)
# size(df)
# describe(df)
# print(df)
delete!(df, [:FILTER,:QUAL, :SOMATIC,:QSI_NT,:TQSI_NT,:SGT])
# head(df)

for c in names(df)
    #println(c)
    k=isna(df[c])
    if mean(k)==1.0
		delete!(df, c)
		continue
	end
	if isa(df[c],DataArray{String,1})
      df[k,c] = ""
	end	
end
#rename!(df, [:normal_tumor_alt_count, :normal_tumor_ref_count], [:t_alt_count, :t_ref_count])
head(df)
# keep only CHRO = [1-Y]
df=df[find(in.(df[:CHRO], [["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]])),:]

# sort by CHRO and position
a= df[:CHRO]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    #a=map(x -> replace(x,r"[MT]", "25"), a)
    a=map(x -> parse(Int32,x), a)
end
df[:a] = a
# print(df)
sort!(df, cols = [:a, :POS])
delete!(df, [:a])


# Strelka 'tier' for QSS -
# hypothesis: TQSS must be strelka's favorite tier for each particular call 
nt=df[:TQSI]
n=length(nt)

q=df[:NORMAL_TIR]
qo=df[:NORMAL_TOR]
NDP=[df[:NORMAL_DP] df[:NORMAL_DP2]]
NI1= fill("0",n)
NO1= fill("0",n)
ND1= fill(0,n)
q = map(x -> split(x,":"), q)
qo = map(x -> split(x,":"), qo)
for i = 1:n
	print(i,"\n")
	NI1[i]=q[i][nt[i]]
	ND1[i]=NDP[i,nt[i]]
	NO1[i]=qo[i][nt[i]]
end
NI1=map(x -> parse(Int32,x), NI1)
NO1=map(x -> parse(Int32,x), NO1)

NR1=ND1-NI1  # NORMAL_TOR can sometimes exceed depth so don't subtract NO1


q=df[:TUMOR_TIR]
qo=df[:TUMOR_TOR]
TDP=[df[:TUMOR_DP] df[:TUMOR_DP2]]
TI1= fill("0",n)
TO1= fill("0",n)
TD1= fill(0,n)
q = map(x -> split(x,":"), q)
qo = map(x -> split(x,":"), qo)
for i = 1:n
	print(i,"\n")
	TI1[i]=q[i][nt[i]]
	TD1[i]=TDP[i,nt[i]]
	TO1[i]=qo[i][nt[i]]
end
TI1=map(x -> parse(Int32,x), TI1)
TO1=map(x -> parse(Int32,x), TO1)

TR1=TD1-TI1  # TUMOR_TOR can sometimes exceed depth so don't subtract NO1


df[:n_ref_count] = NR1
df[:n_alt_count] = NI1
df[:t_ref_count] = TR1
df[:t_alt_count] = TI1
df[:QS] = df[:QSI]


open(file3a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Array,df), '\t')
end

df4=df




df= join(df1, df2, on =  [:CHRO, :POS, :REF, :ALT], kind = :outer)
#df= join(df, df3, on =  [:CHRO, :POS, :REF, :ALT], kind = :outer)
#df= join(df, df4, on =  [:CHRO, :POS, :REF, :ALT], kind = :outer)

# label M1 with M1 flag
df[:M1]=~isna(df[:n_alt_count])
df[:M2]=~isna(df[:TLOD])
#df[:ST]=(df[:QSS]>-1)|(df[:QSI]>-1)

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


df12=df

df= join(df3, df4, on =  [:CHRO, :POS, :REF, :ALT], kind = :outer)

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
df[:QS]=df[:QSS]
df[:TQS]=df[:TQSS]

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
df[k,:QS]=df[k,:QSI]
df[k,:TQS]=df[k,:TQSI]
df[k,:n_alt_count]=df[k,:n_alt_count_1]
df[k,:n_ref_count]=df[k,:n_ref_count_1]
df[k,:t_alt_count]=df[k,:t_alt_count_1]
df[k,:t_ref_count]=df[k,:t_ref_count_1]
df[k,:NORMAL_DP]=df[k,:NORMAL_DP_1]
df[k,:TUMOR_DP]=df[k,:TUMOR_DP_1]

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
df[k,:QS]=df[k,:QSI]
df[k,:TQS]=df[k,:TQSI]
df[k,:n_alt_count]=df[k,:n_alt_count_1]
df[k,:n_ref_count]=df[k,:n_ref_count_1]
df[k,:t_alt_count]=df[k,:t_alt_count_1]
df[k,:t_ref_count]=df[k,:t_ref_count_1]
df[k,:NORMAL_DP]=df[k,:NORMAL_DP_1]
df[k,:TUMOR_DP]=df[k,:TUMOR_DP_1]


# add build column
df[:build]=fill("37",size(df,1))
df[:tumor_barcode]=fill(tumor_id,size(df,1))
df[:normal_barcode]=fill(normal_id,size(df,1))
df[:judgement]=fill("KEEP",size(df,1))

rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])

delete!(df, [:NT,:TQSS,:QSS,:NORMAL_AU,:NORMAL_SDP,:NORMAL_SUBDP,:NORMAL_TU,:NORMAL_GU,:NORMAL_CU,:NORMAL_FDP,:TUMOR_AU ,:TUMOR_SDP,:TUMOR_SUBDP,:TUMOR_TU ,:TUMOR_GU ,:TUMOR_CU ,:TUMOR_FDP])
delete!(df, [:QSI,:IC,:NT_1,:IHP,:TQSI,:RU,:OVERLAP,:RC,:NORMAL_DP50,:NORMAL_DP_1,:NORMAL_FDP50,:NORMAL_SUBDP50,:NORMAL_TOR,:NORMAL_TIR,:NORMAL_DP2,:NORMAL_TAR,:TUMOR_DP50,:TUMOR_DP_1,:TUMOR_FDP50,:TUMOR_SUBDP50,:TUMOR_TOR,:TUMOR_TIR,:TUMOR_DP2,:TUMOR_TAR])
delete!(df, [:n_ref_count_1,:n_alt_count_1,:t_ref_count_1,:t_alt_count_1,:QS_1,:ID,:ID_1])
df34=df

df= join(df12, df34, on =  [:chr, :start, :ref_allele, :alt_allele], kind = :outer)

kST=find(isna(df[:n_alt_count]))
df[kST,:n_alt_count]=df[kST,:n_alt_count_1]
df[kST,:n_ref_count]=df[kST,:n_ref_count_1]
df[kST,:t_alt_count]=df[kST,:t_alt_count_1]
df[kST,:t_ref_count]=df[kST,:t_ref_count_1]
df[kST,:end]=df[kST,:end_1]
df[kST,:tumor_barcode]=df[kST,:tumor_barcode_1]
df[kST,:normal_barcode]=df[kST,:normal_barcode_1]
df[kST,:judgement]=df[kST,:judgement_1]
df[:tumor_f]=df[:t_alt_count]./(df[:t_alt_count]+df[:t_ref_count])

df[:STRLK]=~isna(df[:TQS])

# add build column
df[:build]=fill("37",size(df,1))
df[:tumor_barcode]=fill(tumor_id,size(df,1))
df[:normal_barcode]=fill(normal_id,size(df,1))
df[:judgement]=fill("KEEP",size(df,1))

k=find(isna(df[:NORMAL_DP]))
df[k,:NORMAL_DP]=df[k,:n_alt_count]+df[k,:n_ref_count]
df[k,:TUMOR_DP]=df[k,:t_alt_count]+df[k,:t_ref_count]

k=find(isna(df[:M1]))
df[k,:M1]=false
df[k,:M2]=false


delete!(df, [:n_ref_count_1,:n_alt_count_1,:t_ref_count_1,:t_alt_count_1])
delete!(df, [:end_1,:build_1,:tumor_barcode_1,:normal_barcode_1,:judgement_1])


#from = ["ID_1","DB","HCNT","MAX_ED","RPA","STR","RU","NLOD","TUMOR_ALT_F1R2","MIN_ED","ECNT","TLOD","TUMOR_PID","TUMOR_AF","TUMOR_FOXOG","TUMOR_QSS","TUMOR_ALT_F2R1","TUMOR_AD_REF","TUMOR_AD_ALT","TUMOR_REF_F2R1","TUMOR_REF_F1R2","NORMAL_AF","NORMAL_AD_REF","NORMAL_AD_ALT"]
#println(from)
#
#for c in names(df)
#    c1 = string(c)
#    println(c1)
#    if length(find(Bool[contains(c1,i) for i in from]))>0
#        m2c = Symbol(string("STRLK",c1))
#        rename!(df,c,m2c)
#    end
#end


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
    maf[k,c]=""
    k=find(map(x->x=="NA",maf[c]))
    maf[k,c] = ""
    k=find(map(x->uppercase(x)=="TRUE",maf[c]))
    maf[k,c] = "1"
    k=find(map(x->uppercase(x)=="FALSE",maf[c]))
    maf[k,c] = "0"
end

a= maf[:chr]
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    a=map(x -> parse(Int32,x), a)
end
maf[:a] = a
sort!(maf, cols = [:a, :start])
delete!(maf, [:a])


open(maflite, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
