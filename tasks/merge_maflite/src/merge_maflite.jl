#!/usr/local/bin/julia
# ARGS=["REBC-AC8L-TP","REBC-AC8L-NT","sample.mutect.maflite.txt","REBC-AC8L-TP-NT.m2_maflite.tsv","REBC-AC8L-TP-NT.Strelka_maflite.tsv","REBC-AC8L-TP-NT.SvABA_maflite.tsv","M1","M2","STRELKA","SVABA","REBC-AC8L-TP-NT.merged_maflite.tsv"]
using DataFrames
tumor_id=ARGS[1]
normal_id=ARGS[2]
file1=ARGS[3]
file2=ARGS[4]
file3=ARGS[5]
file4=ARGS[6]
lab1=ARGS[7]
lab2=ARGS[8]
lab3=ARGS[9]
lab4=ARGS[10]
merged_maflite=ARGS[11]

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

maflite_fields=["build","chr","start","end","ref_allele","tum_allele1","tum_allele2","tumor_barcode","normal_barcode","t_alt_count","t_ref_count","judgement","n_alt_count","n_ref_count","tumor_f"]
maflite_symbols=map(x->Symbol(x),maflite_fields)


if isfile(file1)
	df1 = readtable(file1,separator ='\t')
	if Symbol("_end") in names(df1)
		rename!(df1, [:_end], [:end])
	end
	for c in names(df1)
	    if ~isString(df1[1,c])
    	    df1[c] = map(x -> string(x),df1[c])
    	end
    end

    for c in names(df1)
    	c1 = string(c)
    	if length(find(Bool[contains(c1,i) for i in maflite_fields]))<1
        	m2c = Symbol(string(lab1,"_",c1))
        	rename!(df1,c,m2c)
		end
	end
    
    df1[Symbol(lab1)]=fill("1",size(df1[:chr]))
    df=df1
end

if isfile(file2)
	df2 = readtable(file2)
	for c in names(df2)
	    if ~isString(df2[1,c])
    	    df2[c] = map(x -> string(x),df2[c])
    	end
    end
    if Symbol("_end") in names(df2)
		rename!(df2, [:_end], [:end])
	end
	if ! (Symbol("tum_allele1") in names(df2))
		df2[:tum_allele1]=df2[:ref_allele]
		df2[:tum_allele2]=df2[:alt_allele]
	end

	for c in names(df2)
    	c1 = string(c)
    	if length(find(Bool[contains(c1,i) for i in maflite_fields]))<1
        	m2c = Symbol(string(lab2,"_",c1))
        	rename!(df2,c,m2c)
		end
	end
	
    df2[Symbol(lab2)]=fill("1",size(df2[:chr]))
	if isdefined(:df)
	   df[Symbol(lab2)]=fill("0",size(df[:chr]))
	   df=[df; df2]
	else
	   df=df2
	end
end
if isfile(file3)
	df3 = readtable(file3)
	for c in names(df3)
	    if ~isString(df3[1,c])
    	    df3[c] = map(x -> string(x),df3[c])
    	end
    end
	if Symbol("_end") in names(df3)
		rename!(df3, [:_end], [:end])
	end
	if ! (Symbol("tum_allele1") in names(df3))
		df3[:tum_allele1]=df3[:ref_allele]
		df3[:tum_allele2]=df3[:alt_allele]
	end

	for c in names(df3)
    	c1 = string(c)
    	if length(find(Bool[contains(c1,i) for i in maflite_fields]))<1
        	m2c = Symbol(string(lab3,"_",c1))
        	rename!(df3,c,m2c)
		end
	end
	
    df3[Symbol(lab3)]=fill("1",size(df3[:chr]))
	if isdefined(:df)
	   df[Symbol(lab3)]=fill("0",size(df[:chr]))
	   df=[df; df3]
	else
	   df=df3
	end
end
if isfile(file4)
	df4 = readtable(file4)
	for c in names(df4)
	    if ~isString(df4[1,c])
    	    df4[c] = map(x -> string(x),df4[c])
    	end
    end
	if Symbol("_end") in names(df4)
		rename!(df4, [:_end], [:end])
	end
	if ! (Symbol("tum_allele1") in names(df4))
		df4[:tum_allele1]=df4[:ref_allele]
		df4[:tum_allele2]=df4[:alt_allele]
	end

	for c in names(df4)
    	c1 = string(c)
    	if length(find(Bool[contains(c1,i) for i in maflite_fields]))<1
        	m2c = Symbol(string(lab4,"_",c1))
        	rename!(df4,c,m2c)
		end
	end
	
    df4[Symbol(lab4)]=fill("1",size(df4[:chr]))
	if isdefined(:df)
	   df[Symbol(lab4)]=fill("0",size(df[:chr]))
	   df=[df; df4]
	else
	   df=df4
	end
end

# fix fields 

df[:build]=fill("37",size(df,1))
df[:tumor_barcode]=fill(tumor_id,size(df,1))
df[:normal_barcode]=fill(normal_id,size(df,1))
df[:judgement]=fill("KEEP",size(df,1))

t_alt_count = map(x -> parse(Float64,x),df[:t_alt_count])
t_ref_count = map(x -> parse(Float64,x),df[:t_ref_count])
df[:tumor_f]= map(x->@sprintf("%.4f",x),t_alt_count./(t_alt_count+t_ref_count))

# sort 
a = df[:chr]
p1 = map(x -> parse(Int64,x), df[:start])
p2 = map(x -> parse(Int64,x), df[:end])
if !isa(a[1],Int)
    a=map(x -> replace(x,r"[X]", "23"), a)
    a=map(x -> replace(x,r"[Y]", "24"), a)
    a=map(x -> parse(Int64,x), a)
end
df[:a] = a
df[:p1] = p1
df[:p2] = p2
#print(df)
sort!(df, cols = [:a, :p1])

# convert NA's to ""
for c in names(df)
    k=find(map(x -> isna(x),df[c]))
    df[k,c] = ""
    if c in [Symbol(lab1), Symbol(lab2), Symbol(lab3), Symbol(lab4)]
	   k=find(map(x -> x=="",df[c]))
 	   df[k,c] = "0"
    end    
end	

# cluster
x1=df[:a]*1000000000+df[:p1]
x2=df[:a]*1000000000+df[:p2]
kdel = find(map(x-> x=="-", df[:ref_allele]))
x1[kdel]=x1[kdel]-1
x2[kdel]=x2[kdel]+1
kins = find(map(x-> x=="-", df[:tum_allele2]))
x1[kins]=x1[kins]-1
x2[kins]=x2[kins]+1

n=length(df[:a])
g=1
df[:cluster]=fill(g,size(df[:chr]))
for i = 2:n
	if (x1[i]>x2[i-1])
	  g=g+1
	end
    df[i,:cluster]=g
end

# q=countmap(df[:cluster]) ; q1=map(x->Int32(x),values(q)); countmap(q1)

fields=names(df)
sfields=map(x -> string(x),names(df))

labs=[lab1,lab2,lab3,lab4]
slabs=map(x -> Symbol(x),labs)
			
df[:NALG]=fill(1,size(df[:chr]))

for i = 1:g
	
    k=find(map(x-> x==i,df[:cluster]))
    #println(k)
    if length(k)>1
    	if (any(map(x->x=="-",df[k,:ref_allele])) | any(map(x->x=="-",df[k,:tum_allele2]) ) )
	  		println(i)
      		println(df[k,[:chr,:start,:ref_allele,:tum_allele2,:t_alt_count,:tumor_f,:M1,:M2,:STRELKA,:SVABA]])
    	end
    	
    	# algorithm pecking order as input lab1-lab4 
    	# highest t_alt_count wins 
    	t_alt_count = map(x->parse(Int32,x),df[k,:t_alt_count])
    	k1=k[1]
    	kx=k[2:end]
    	kk=k[findmax(t_alt_count)[2]]
    	
    	# each algorithm can appear only once
    	#df[k,slabs]
    	for a=1:4
    	   if (sum(map(x -> x=="1",df[k,slabs[a]])))>1
    	   	 
    	   	 println(i)
        	 println(df[k,[:chr,:start,:ref_allele,:tum_allele2,:t_alt_count,:tumor_f,:M1,:M2,:STRELKA,:SVABA]])

    	     ka=kx[find(map(x -> x=="1",df[kx,slabs[a]]))][1]
    	     df[ka,:judgement]="ALG_REDUNDANT" 
    	     # remove ka from k, redefine t_alt_count, k1,kx, kk 
    	     filter!(e->e≠ka,k)
    	     if length(k)<2
    	        continue
    	     end
    	     t_alt_count = map(x->parse(Int32,x),df[k,:t_alt_count])
    	     k1=k[1]
    	     kx=k[2:end]
	         kk=k[findmax(t_alt_count)[2]]    	      
	      end  
    	end
    	if length(k)<2
            continue
    	end
    	    	
    	df[k1,maflite_symbols]=df[kk,maflite_symbols]
        df[k1,:NALG]=length(k)
        df[kx,:judgement]="REDUNDANT" 
        for lab=labs
			slab=Symbol(lab)
	    	#println(slab)
	    	mfields=fields[find(map(x -> searchindex(x,lab)==1,sfields))]
	    	kalg=k[find(map(x -> x=="1",df[k,slab]))] 
	    	if ~isempty(kalg)	    	    
	    		df[k1,mfields]=df[kalg,mfields]
			end
	    	
    	end
    	
    end
end

k=find(map(x -> x=="KEEP",df[:judgement]))

delete!(df, [:a,:p1,:p2,:cluster])

maf=df[k,:]
for c in names(maf)
    #println(c)
    if ~isString(maf[1,c])
        maf[c] = map(x -> string(x),maf[c])
    end
end





open(merged_maflite, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
