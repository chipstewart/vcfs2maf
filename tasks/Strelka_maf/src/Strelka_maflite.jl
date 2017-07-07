 #!/usr/local/bin/julia
#ARGS
# tumor_id="T"; normal_id="N"; file1="tmp1.tsv"; file2="tmp2.tsv"; maflite="maflite.tsv"
# ARGS=["THCA-BJ-A28T-TP","THCA-BJ-A28T-NB","THCA-BJ-A28T-TP-NB.StrelkaSNV.tsv","THCA-BJ-A28T-TP-NB.StrelkaINDEL.tsv","THCA-BJ-A28T-TP-NB.Strelka.maflite.tsv"]
using DataFrames

file1a="vcf_1.tsv"
file2a="vcf_2.tsv"

tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
file2=ARGS[4]
maflite=ARGS[5]


#Strelka SNV vcf
df = readtable(file1)
# size(df)
# describe(df)
# print(df)
delete!(df, [:FILTER,:QUAL,:ID,:NT,:SOMATIC,:QSS_NT,:TQSS_NT,:SGT])
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
# if CHRO is numberical -> convert to string
if (typeof(df[:CHRO])==typeof(DataArray(fill(34,2))))
	 df[:CHRO]=map(x -> string(x),df[:CHRO])
end
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

q=df[:TUMOR_GU]
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


open(file1a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Array,df), '\t')
end

df1=df

#Strelka INDEL vcf
df = readtable(file2)
# size(df)
# describe(df)
# print(df)
delete!(df, [:FILTER,:QUAL,:ID,:NT,:SOMATIC,:QSI_NT,:TQSI_NT,:SGT,:SVTYPE,:OVERLAP])
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
if (typeof(df[:CHRO])==typeof(DataArray(fill(34,2))))
	 df[:CHRO]=map(x -> string(x),df[:CHRO])
end
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


open(file2a, "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Array,df), '\t')
end

df2=df





df= join(df1, df2, on =  [:CHRO, :POS, :REF, :ALT], kind = :outer)

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

rename!(df, [:CHRO,:POS, :REF, :ALT], [:chr, :start, :ref_allele,:alt_allele])

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

n=length(maf[:chr])
maf[:NORMAL_ACGT_TIR_TOR]=fill("",n)
maf[:TUMOR_ACGT_TIR_TOR]=fill("",n)
for i = 1:n
	print(i,"\n")
	maf[:NORMAL_ACGT_TIR_TOR][i]=string(maf[:NORMAL_AU][i],';',df[:NORMAL_CU][i],';',maf[:NORMAL_GU][i],';',df[:NORMAL_TU][i],';',df[:NORMAL_TIR][i],';',df[:NORMAL_TOR][i])
	maf[:TUMOR_ACGT_TIR_TOR][i]=string(maf[:TUMOR_AU][i],';',df[:TUMOR_CU][i],';',maf[:TUMOR_GU][i],';',df[:TUMOR_TU][i],';',df[:TUMOR_TIR][i],';',df[:TUMOR_TOR][i])
end

delete!(maf, [:NORMAL_AU,:NORMAL_SDP,:NORMAL_SUBDP,:NORMAL_TU,:NORMAL_GU,:NORMAL_CU,:NORMAL_DP_1])
delete!(maf, [:TUMOR_AU,:TUMOR_SDP,:TUMOR_SUBDP,:TUMOR_TU,:TUMOR_GU,:TUMOR_CU,:TUMOR_DP_1])
delete!(maf, [:TQSS,:QSS,:TQSI,:QSI,:IC,:IHP])
delete!(maf, [:NORMAL_DP50,:NORMAL_FDP50,:NORMAL_SUBDP50,:NORMAL_TOR,:NORMAL_TIR,:NORMAL_DP2,:NORMAL_TAR])
delete!(maf, [:TUMOR_DP50,:TUMOR_FDP50,:TUMOR_SUBDP50,:TUMOR_TOR,:TUMOR_TIR,:TUMOR_DP2,:TUMOR_TAR])
delete!(maf, [:n_ref_count_1,:n_alt_count_1,:t_ref_count_1,:t_alt_count_1,:QS_1])

rename!(df, [:NORMAL_DP,:NORMAL_FDP,:TUMOR_DP, :TUMOR_FDP], [:STRELKA_NORMAL_DP, :STRELKA_NORMAL_FDP, :STRELKA_TUMOR_DP, :STRELKA_TUMOR_FDP])
rename!(df, [:QS,:RU,:RC, :TQS], [:STRELKA_QS, :STRELKA_RU, :STRELKA_RC, :STRELKA_TQS])
rename!(df, [:NORMAL_ACGT_TIR_TOR,:TUMOR_ACGT_TIR_TOR], [:STRELKA_NORMAL_ACGT_TIR_TOR, :STRELKA_TUMOR_ACGT_TIR_TOR])




open(maflite, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
