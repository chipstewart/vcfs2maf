#!/usr/local/bin/julia
# ARGS=["REBC-AC8L-TP","REBC-AC8L-NT","sample.mutect.maflite.txt","REBC-AC8L-TP-NT.m2_maflite.tsv","REBC-AC8L-TP-NT.Strelka_maflite.tsv","REBC-AC8L-TP-NT.SvABA_maflite.tsv","REBC-AC8L-TP-NT.snowman_maflite.tsv","M1","M2","STRELKA","SVABA","Snowman",REBC-AC8L-TP-NT.merged_maflite.tsv"]
# ARGS=["REBC-AC8R-TP","REBC-AC8R-NB","REBC-AC8R-TP-NB.m1_maflite.tsv","REBC-AC8R-TP-NB.m2_maflite.tsv","REBC-AC8R-TP-NB.Strelka_maflite.tsv","REBC-AC8R-TP-NB.SvABA_maflite.tsv","REBC-AC8R-TP-NB.snowman_maflite.tsv","M1","M2","STRELKA","SVABA","Snowman","REBC-AC8R-TP-NB.merged_maflite.tsv"]
# ARGS=["REBC-AF7Y-TTP1-A-1-1-D-A649-36","SC208303","REBC-AF7Y-TP-NB.m1_maflite.tsv","REBC-AF7Y-TP-NB.m2_maflite.tsv","REBC-AF7Y-TP-NB.Strelka_maflite.tsv","REBC-AF7Y-TP-NB.SvABA_maflite.tsv","REBC-AF7Y-TP-NB.snowman_maflite.tsv","REBC-AF7Y-TP-NB.Strelka2_maflite.tsv","M1","M2","Strelka1","SVABA","Snowman","Strelka2","REBC-AF7Y-TP-NB.merged_maflite.tsv"]
# ARGS=["REBC-ACBQ-TTM1-A-1-1-D-A76R-36","REBC-ACBQ-NT1-A-1-1-D-A49W-36","REBC-ACBQ-TM-NT.m1_maflite.tsv","REBC-ACBQ-TM-NT.m2_maflite.tsv","REBC-ACBQ-TM-NT.Strelka_maflite.tsv","REBC-ACBQ-TM-NT.Strelka2_maflite.tsv","REBC-ACBQ-TM-NT.SvABA_maflite.tsv","REBC-ACBQ-TM-NT.snowman_maflite.tsv","M1","M2","Strelka","Strelka2","SvABA","Snowman","REBC-ACBQ-TM-NT.merged.maflite.tsv"]


using DataFrames
using Missings
using CSV
using Printf
using DelimitedFiles


print(ARGS)

tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
file2=ARGS[4]
file3=ARGS[5]
file4=ARGS[6]
file5=ARGS[7]
file6=ARGS[8]
lab1=ARGS[9]
lab2=ARGS[10]
lab3=ARGS[11]
lab4=ARGS[12]
lab5=ARGS[13]
lab6=ARGS[14]
merged_maflite=ARGS[15]

isString(x::Number)=false
isString(x::Missing)=false
isString(x::AbstractString)=true

maflite_fields=["build","chr","start","end","ref_allele","alt_allele","tumor_barcode","normal_barcode","t_alt_count","t_ref_count","judgement","n_alt_count","n_ref_count","tumor_f"]
maflite_symbols=map(x->Symbol(x),maflite_fields)


if isfile(file1)&&(lab1!="-")
    #df1 = readtable(file1,separator ='\t')
    df1 = CSV.File(file1,delim ='\t') |> DataFrame
    if Symbol("_end") in names(df1)
        rename!(df1, [:_end], [:end])
    end
    for c in names(df1)
        if ~isString(df1[1,c])
            df1[!,c] = map(x -> string(x),df1[!,c])
        end
    end

    for c in names(df1)
        c1 = string(c)
        if length(findall(Bool[occursin(c1,i) for i in maflite_fields]))<1
            m2c = Symbol(string(lab1,"_",c1))
            rename!(df1,c => m2c)
        end
    end
    
    df1[!,Symbol(lab1)]=fill("1",size(df1[!,:chr]))
    df=df1
end

show(df)

if isfile(file2)&&(lab2!="-")
    df2 =  CSV.File(file2,delim ='\t') |> DataFrame
    for c in names(df2)
        if ~isString(df2[1,c])
            df2[!,c] = map(x -> string(x),df2[!,c])
        end
    end
    if Symbol("_end") in names(df2)
        rename!(df2, [:_end] => [:end])
    end

    for c in names(df2)
        c1 = string(c)
        if length(findall(Bool[occursin(c1,i) for i in maflite_fields]))<1
            m2c = Symbol(string(lab2,"_",c1))
            rename!(df2,c => m2c)
        end
    end
    
    df2[!,Symbol(lab2)]=fill("1",size(df2[!,:chr]))
    if @isdefined df
       df[!,Symbol(lab2)]=fill("0",size(df[!,:chr]))
       #df=[df; df2]
       df=outerjoin(df, df2, on = intersect(names(df), names(df2)))
    else
       df=df2
    end
end

if isfile(file3)&&(lab3!="-")
    df3 = CSV.File(file3,delim ='\t') |> DataFrame
    for c in names(df3)
        if ~isString(df3[1,c])
            df3[!,c] = map(x -> string(x),df3[!,c])
        end
    end
    if Symbol("_end") in names(df3)
        rename!(df3, [:_end] => [:end])
    end

    for c in names(df3)
        c1 = string(c)
        if length(findall(Bool[occursin(c1,i) for i in maflite_fields]))<1
            m2c = Symbol(string(lab3,"_",c1))
            rename!(df3,c => m2c)
        end
    end
    
    df3[!,Symbol(lab3)]=fill("1",size(df3[!,:chr]))
    if @isdefined df 
       df[!,Symbol(lab3)]=fill("0",size(df[!,:chr]))
       df=outerjoin(df, df3, on = intersect(names(df), names(df3)))   
    else
       df=df3
    end
end

if isfile(file4)&&(lab4!="-")
    df4 = CSV.File(file4,delim ='\t') |> DataFrame
    for c in names(df4)
        if ~isString(df4[1,c])
            df4[!,c] = map(x -> string(x),df4[!,c])
        end
    end
    if Symbol("_end") in names(df4)
        rename!(df4, [:_end] => [:end])
    end

    for c in names(df4)
        c1 = string(c)
        if length(findall(Bool[occursin(c1,i) for i in maflite_fields]))<1
            m2c = Symbol(string(lab4,"_",c1))
            rename!(df4,c => m2c)
        end
    end
    
    df4[!,Symbol(lab4)]=fill("1",size(df4[!,:chr]))
    if @isdefined df
       df[!,Symbol(lab4)]=fill("0",size(df[!,:chr]))
       df=outerjoin(df, df4, on = intersect(names(df), names(df4)))  
    else
       df=df4
    end
end

if isfile(file5)&&(lab5!="-")
    df5 = CSV.File(file5,delim ='\t') |> DataFrame
    for c in names(df5)
        if ~isString(df5[1,c])
            df5[!,c] = map(x -> string(x),df5[!,c])
        end
    end
    if Symbol("_end") in names(df5)
        rename!(df5, [:_end] => [:end])
    end

    for c in names(df5)
        c1 = string(c)
        if length(findall(Bool[occursin(c1,i) for i in maflite_fields]))<1
            m2c = Symbol(string(lab5,"_",c1))
            rename!(df5,c => m2c)
        end
    end
    
    df5[!,Symbol(lab5)]=fill("1",size(df5[!,:chr]))
    if @isdefined df
       df[!,Symbol(lab5)]=fill("0",size(df[!,:chr]))
       df=outerjoin(df, df5, on = intersect(names(df), names(df5)))  
    else
       df=df5
    end
end

if isfile(file6)&&(lab6!="-")
    df6 = CSV.File(file6,delim ='\t') |> DataFrame
    for c in names(df6)
        if ~isString(df6[1,c])
            df6[!,c] = map(x -> string(x),df6[!,c])
        end
    end
    if Symbol("_end") in names(df6)
        rename!(df6, [:_end] => [:end])
    end

    for c in names(df6)
        c1 = string(c)
        if length(findall(Bool[occursin(c1,i) for i in maflite_fields]))<1
            m2c = Symbol(string(lab6,"_",c1))
            rename!(df6,c => m2c)
        end
    end
    
    df6[!,Symbol(lab6)]=fill("1",size(df6[!,:chr]))
    if @isdefined df
       df[!,Symbol(lab6)]=fill("0",size(df[!,:chr]))
       df=outerjoin(df, df6, on = intersect(names(df), names(df6))) 
    else
       df=df6
    end
end

show(df)

# fix fields 
k=findall(map(x -> x!="MT",df[!,:chr]))
df=df[k,:]
k=findall(map(x -> x!="M",df[!,:chr]))
df=df[k,:]

df[!,:build]=fill("37",size(df,1))
df[!,:tumor_barcode]=fill(tumor_id,size(df,1))
df[!,:normal_barcode]=fill(normal_id,size(df,1))
df[!,:judgement]=fill("KEEP",size(df,1))

t_alt_count = map(x -> parse(Float64,x),df[!,:t_alt_count])
t_ref_count = map(x -> parse(Float64,x),df[!,:t_ref_count])
df[!,:tumor_f]= map(x->@sprintf("%.4f",x),t_alt_count./(t_alt_count+t_ref_count))

# sort 
a = df[!,:chr]
p1 = map(x -> parse(Int64,x), df[!,:start])
p2 = map(x -> parse(Int64,x), df[!,:end])
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]" => "23"), a)
    a=map(x -> replace(x,r"[Y]" => "24"), a)
    a=map(x -> parse(Int64,x), a)
end
df[!,:a] = a
df[!,:p1] = p1
df[!,:p2] = p2
#print(df)
sort!(df, [:a, :p1])

# convert NA's to ""
#df0=deepcopy(df)
for c1 in names(df)
    c=Symbol(c1)
    println(c)
    if c in [Symbol(lab1), Symbol(lab2), Symbol(lab3), Symbol(lab4),Symbol(lab5),Symbol(lab6)]
       println(c)
       replace!(df[!,c], missing => "0")
    end    
    replace!(df[!,c], missing => "")
end 

# cluster
x1=df[!,:a]*1000000000+df[!,:p1]
x2=df[!,:a]*1000000000+df[!,:p2]
kdel = findall(map(x-> x=="-", df[!,:ref_allele]))
x1[kdel]=x1[kdel].-1
x2[kdel]=x2[kdel].+1
kins = findall(map(x-> x=="-", df[!,:alt_allele]))
x1[kins]=x1[kins].-1
x2[kins]=x2[kins].+1

n=length(df[!,:a])
g=1
df[!,:cluster]=fill(g,size(df[!,:chr]))
for i = 2:n
    if (x1[i]>x2[i-1])
      global g=g+1
    end
    df[i,:cluster]=g
end

# q=countmap(df[:cluster]) ; q1=map(x->Int32(x),values(q)); countmap(q1)

fields=names(df)
sfields=map(x -> string(x),names(df))
print(sfields)

labs=[lab1,lab2,lab3,lab4,lab5,lab6]
#deleteat!(labs, findin(labs, ["-"]))
labs=filter!(x->x≠"-",labs)
NA=length(labs)
slabs=map(x -> Symbol(x),labs)
            
df[!,:NALG]=fill(0,size(df[!,:chr]))


fg=[:chr,:start,:t_alt_count,:t_ref_count,slabs[1],slabs[2],slabs[3],slabs[4],slabs[5],slabs[6],:judgement]

for g1 = 1:g
    kg=findall(map(x-> x==g1,df[!,:cluster]))
    #println(k)
    if length(kg)<2
        continue
    end
    println("g= $g1")
 
    t_alt_count = map(x->parse(Int32,x),df[kg,:t_alt_count])
    k1=kg[findmax(t_alt_count)[2]]
    kx=kg[findall(map(x->x≠k1,kg))]
    for kx1 in kx
         slab1=findall(map(x -> x=="1",df[kx1,slabs]))[1]
         println(slab1)
         df[k1,slab1]="1"
         df[kx1,:judgement]="ALG_REDUNDANT" 

         pre1=String(slab1) *"_"
         salg=sfields[findall( x -> occursin(pre1, x), sfields)]

         for salg1 in salg
             df[k1,Symbol(salg1)]=df[kx1,Symbol(salg1)]
         end
    end
    show(df[kg,fg])
end


df[!,:NALG]=sum(eachcol(parse.(Int64, df[:, slabs])))



for c in names(df)
    if ~isString(df[1,c])
        df[!,c] = map(x -> string(x),df[!,c])
    end
end

k=findall(map(x -> x=="KEEP",df[!,:judgement]))
maf=df[k,:]
select!(maf, Not([:a,:p1,:p2,:cluster]))


open("premerged_maflite.tsv", "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Matrix,df), '\t')
end

open(merged_maflite, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Matrix,maf), '\t')
end


