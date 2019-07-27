#!/usr/local/bin/julia
# ARGS=["REBC-AC8L-TP","REBC-AC8L-NT","sample.mutect.maflite.txt","REBC-AC8L-TP-NT.m2_maflite.tsv","REBC-AC8L-TP-NT.Strelka_maflite.tsv","REBC-AC8L-TP-NT.SvABA_maflite.tsv","REBC-AC8L-TP-NT.snowman_maflite.tsv","M1","M2","STRELKA","SVABA","Snowman",REBC-AC8L-TP-NT.merged_maflite.tsv"]
#ARGS=["REBC-AC9F-NT1-A-1-1-D-A49U-36","REBC-AC9F-TTP1-A-1-1-D-A49U-36","/Users/stewart/Downloads/REBC-AC9F-NT-TP.m1_maflite.tsv","/Users/stewart/Downloads/REBC-AC9F-NT-TP.m2.deTiN.maflite.tsv","/Users/stewart/Downloads/GDAC_FC_NULL","/Users/stewart/Downloads/REBC-AC9F-NT-TP.Strelka2_maflite.tsv","/Users/stewart/Downloads/REBC-AC9F-NT-TP.SvABA_maflite.tsv","/Users/stewart/Downloads/GDAC_FC_NULL","M1","M2","-","Strelka2","SvABA","-","REBC-AC9F-NT-TP.merged.maflite.tsv","37"]
#ARGS=["REBC-AF65-NT1-A-1-1-D-A649-36","REBC-AF65-TTP1-A-1-1-D-A649-36","/Users/stewart/Downloads/REBC-AF65-NT-TP.m1_maflite.tsv","/Users/stewart/Downloads/REBC-AF65-NT-TP.m2.deTiN.maflite.tsv","/Users/stewart/Downloads/GDAC_FC_NULL","/Users/stewart/Downloads/REBC-AF65-NT-TP.Strelka2_maflite.tsv","/Users/stewart/Downloads/REBC-AF65-NT-TP.SvABA_maflite.tsv","/Users/stewart/Downloads/GDAC_FC_NULL","M1","M2","-","Strelka2","SvABA","-","REBC-AF65-NT-TP.merged.maflite.tsv","37"]

using DataFrames
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
build=ARGS[16]

print(ARGS)

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

#maflite_fields=["build","chr","start","end","ref_allele","tum_allele1","tum_allele2","tumor_barcode","normal_barcode","t_alt_count","t_ref_count","judgement","n_alt_count","n_ref_count","tumor_f"]
maflite_fields=["build","chr","start","end","ref_allele","alt_allele","tumor_barcode","normal_barcode","t_alt_count","t_ref_count","judgement","n_alt_count","n_ref_count","tumor_f"]
maflite_symbols=map(x->Symbol(x),maflite_fields)


if isfile(file1)&&(lab1!="-")
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

if isfile(file2)&&(lab2!="-")
	df2 = readtable(file2)
	for c in names(df2)
	    if ~isString(df2[1,c])
    	    df2[c] = map(x -> string(x),df2[c])
    	end
    end
    if Symbol("_end") in names(df2)
		rename!(df2, [:_end], [:end])
	end
#   if ! (Symbol("tum_allele1") in names(df2))
#		df2[:tum_allele1]=df2[:ref_allele]
#		df2[:tum_allele2]=df2[:alt_allele]
#	end

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

if isfile(file3)&&(lab3!="-")
	df3 = readtable(file3)
	for c in names(df3)
	    if ~isString(df3[1,c])
    	    df3[c] = map(x -> string(x),df3[c])
    	end
    end
	if Symbol("_end") in names(df3)
		rename!(df3, [:_end], [:end])
	end
#	if ! (Symbol("tum_allele1") in names(df3))
#		df3[:tum_allele1]=df3[:ref_allele]
#		df3[:tum_allele2]=df3[:alt_allele]
#	end

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

if isfile(file4)&&(lab4!="-")
	df4 = readtable(file4)
	for c in names(df4)
	    if ~isString(df4[1,c])
    	    df4[c] = map(x -> string(x),df4[c])
    	end
    end
	if Symbol("_end") in names(df4)
		rename!(df4, [:_end], [:end])
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


if isfile(file5)&&(lab5!="-")
    df5 = readtable(file5)
    for c in names(df5)
        if ~isString(df5[1,c])
            df5[c] = map(x -> string(x),df5[c])
        end
    end
    if Symbol("_end") in names(df5)
        rename!(df5, [:_end], [:end])
    end

    for c in names(df5)
        c1 = string(c)
        if length(find(Bool[contains(c1,i) for i in maflite_fields]))<1
            m2c = Symbol(string(lab5,"_",c1))
            rename!(df5,c,m2c)
        end
    end
    
    df5[Symbol(lab5)]=fill("1",size(df5[:chr]))
    if isdefined(:df)
       df[Symbol(lab5)]=fill("0",size(df[:chr]))
       df=[df; df5]
    else
       df=df5
    end
end

if isfile(file6)&&(lab6!="-")
    df6 = readtable(file6)
    for c in names(df6)
        if ~isString(df6[1,c])
            df6[c] = map(x -> string(x),df6[c])
        end
    end
    if Symbol("_end") in names(df6)
        rename!(df6, [:_end], [:end])
    end

    for c in names(df6)
        c1 = string(c)
        if length(find(Bool[contains(c1,i) for i in maflite_fields]))<1
            m2c = Symbol(string(lab6,"_",c1))
            rename!(df6,c,m2c)
        end
    end
    
    df6[Symbol(lab6)]=fill("1",size(df6[:chr]))
    if isdefined(:df)
       df[Symbol(lab6)]=fill("0",size(df[:chr]))
       df=[df; df6]
    else
       df=df6
    end
end

# fix fields 
k=find(map(x -> x!="MT",df[:chr]))
df=df[k,:]
k=find(map(x -> x!="M",df[:chr]))
df=df[k,:]

# build from args 
df[:build]=fill(build,size(df,1)) 
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
    if c in [Symbol(lab1), Symbol(lab2), Symbol(lab3), Symbol(lab4),Symbol(lab5),Symbol(lab6)]
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
kins = find(map(x-> x=="-", df[:alt_allele]))
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
print(sfields)

labs=[lab1,lab2,lab3,lab4,lab5,lab6]
deleteat!(labs, findin(labs, ["-"]))
NA=length(labs)
slabs=map(x -> Symbol(x),labs)
			
df[:NALG]=fill(1,size(df[:chr]))

kpre=0
for i = 1:g
    #k=find(map(x-> x==i,df[:cluster])) # very slow when length(df)>500k - but :cluster is monotonic so it's faster to loop
	k=Int64[]
	i1=kpre+1
	while df[i1,:cluster]==i
		k=vcat(k,i1)
		i1=i1+1
		if (i1>n)
		   break
		end
	end	
    kpre=k[end]
    println(i,"/",g,":",k)
    if length(k)>1
    	# if (any(map(x->x=="-",df[k,:ref_allele])) | any(map(x->x=="-",df[k,:tum_allele2]) ) )
	  		# println(i)
            # println(df[k,[:chr,:start,:ref_allele,:tum_allele2,:t_alt_count,:tumor_f,:M1,:M2,:STRELKA,:SVABA]])
    	# end
    	
    	# algorithm pecking order as input lab1-lab4 
    	# highest t_alt_count wins 
    	t_alt_count = map(x->parse(Int32,x),df[k,:t_alt_count])
    	k1=k[1]
    	kx=k[2:end]
    	kk=k[findmax(t_alt_count)[2]]
    	
    	# each algorithm can appear only once
    	#df[k,slabs]
    	for a=1:NA
    	   if (sum(map(x -> x=="1",df[k,slabs[a]])))>1
    	   	 
    	   	 #println(i)
        	 #println(df[k,[:chr,:start,:ref_allele,:tum_allele2,:t_alt_count,:tumor_f,:M1,:M2,:STRELKA,:SVABA]])

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
    if length(k)<1
      println(i,"/",g,":",k)
      error("skipped group in loop over groups")
    end
end


for c in names(df)
    if ~isString(df[1,c])
        df[c] = map(x -> string(x),df[c])
    end
end

k=find(map(x -> x=="KEEP",df[:judgement]))
maf=df[k,:]
delete!(maf, [:a,:p1,:p2,:cluster])

open("premerged_maflite.tsv", "w") do f
    writedlm(f, reshape(names(df), 1, length(names(df))), '\t')
    writedlm(f, convert(Array,df), '\t')
end

open(merged_maflite, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end
