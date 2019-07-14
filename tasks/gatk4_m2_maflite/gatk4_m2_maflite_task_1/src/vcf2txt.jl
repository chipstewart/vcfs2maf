#!/usr/local/bin/julia
# export JULIA_LOAD_PATH="/path/to/module/storage/folder"
#include("/Users/stewart/Projects/jlVCF/src/jlVCF.jl")
# ARGS=["/Users/stewart/Downloads/753_pair.MuTect2.call_stats.vcf","/Users/stewart/PycharmProjects/vcf2txt/out/test.tsv"]

include("jlVCF.jl")
include("VCFIterator.jl")
include("Variant.jl")

narg = 0
f=stdout
max_events=2^31
vcf_file1=ARGS[1]


narg=length(ARGS)
if (narg>1)
	a=ARGS[2]
	#if isnull(tryparse(Int32,a))
	if tryparse(Int32,a)==nothing
		println("vcf:\t",ARGS[1])
		println("tsv:\t",a)
		f = open(a,"w")
	else
	  max_events=parse(Int32,a)
	end
end

delim="\t"
c=fill("", 500)
num=Array{AbstractString,500}
if (max_events<1)
	delim="\n"
	c=Array{AbstractString,500}
	for i in 1:500
	  c[i]=string(i,".\t")
	end
end
#test loading
vc = VCFIterator(vcf_file1) 

print(f,c[1],"CHRO",delim)
print(f,c[2],"POS",delim)
print(f,c[3],"ID",delim)
print(f,c[4],"REF",delim)
print(f,c[5],"ALT",delim)
print(f,c[6],"QUAL",delim)
print(f,c[8],"FILTER",delim)
global n1=8
for x in keys(vc.infoTypes)
    global n1+=1
	print(f,c[n1],x,delim)
end

NS=length(vc.samples)
for s in vc.samples
	for x in keys(vc.formatTypes)    
		if vc.formatTypes[x].Number=="R"
		    global n1+=1
	  		print(f,c[n1],s,"_",x,"_REF",delim)
		    global n1+=1
	  		print(f,c[n1],s,"_",x,"_ALT",delim)
		else
		    global n1+=1
	  		print(f,c[n1],s,"_",x,delim)
		end	
	end
end
print(f,"\n" )

 

count = 0
v = "dmy"
while !eof(vc)
    v = next1(vc)
    if isa(v,Bool)
    	continue
    end
    if sizeof(v)<2
    	continue
    end
    global count += 1
    if count>max_events
	  break
	end
	#print(f,count,"\t")
	print(f,v.CHROM,"\t")
	print(f,v.POS,"\t")
	print(f,join(v.ID,':'),"\t")
	print(f,v.REF,"\t")
	print(f,join(v.ALT,':'),"\t")
	print(f,v.QUAL,"\t")
	print(f,join(v.FILTER,':'),"\t")
	for x in keys(vc.infoTypes)
	  print(f,join((get(v.INFO,x,"")),":"),"\t")
	end	
	for s in vc.samples
		for x in keys(vc.formatTypes)
		  if vc.formatTypes[x].Number=="R"
			print(f,v.FORMAT[s][x][1],"\t")
	  		print(f,v.FORMAT[s][x][2],"\t")
	 	  else
	  		q=get(v.FORMAT[s],x,"")
	  		print(f,join(q,":"),"\t")
	  	  end
	  	end
	end	
	print(f,"\n" )
end

close(vc)
