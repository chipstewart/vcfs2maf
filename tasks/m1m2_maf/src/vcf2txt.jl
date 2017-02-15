#!/usr/local/bin/julia
# export JULIA_LOAD_PATH="/path/to/module/storage/folder"
#include("/Users/stewart/Projects/jlVCF/src/jlVCF.jl")

using jlVCF

narg = 0
f=STDOUT
max_events=2^31
#vcf_file1="/Volumes/64GB_2015/Work/REBC/sample.broad-mutect.DATECODE.somatic.snv_mnv.vcf"
vcf_file1=ARGS[1]

# is this a gz file? 
# f=splitext(vcf_file1)
# if (uppercase(f[2])==".GZ")
#    f0=vcf_file1
#    f1=splitdir(f[1])
#    f1=f[1]
#    f2=splitdir(f1)
#    f2=f2[2]
#    ;zcat $f0 > $f2
#    #run(`zcat' $f0  |> $f2')
#     #run(pipeline(`zcat $f0`, stdout=$f2, stderr="zcat.stderr.log"))
#    vcf_file1=f2
# end

narg=length(ARGS)
#println(narg)
if (narg>1)
	a=ARGS[2]
	if  (narg>1) 
		if isnull(tryparse(Int32,a))
	     	println("vcf:\t",ARGS[1])
	     	println("tsv:\t",a)
	     	f = open(a,"w")
	    else
	      max_events=parse(Int32,a)
	    end
	end
end

delim="\t"
c=fill("", 500)
num=Array{AbstractString}(500)
if (max_events<1)
	delim="\n"
	c=Array{AbstractString}(500)
	for i in 1:500
	  c[i]=dec(i)".\t"
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
n=8
for x in keys(vc.infoTypes)
    n=n+1
	print(f,c[n],x,delim)
end

NS=length(vc.samples)
for s in vc.samples
	for x in keys(vc.formatTypes)    
		if vc.formatTypes[x].Number=="R"
		    n=n+1
	  		print(f,c[n],s,"_",x,"_REF",delim)
		    n=n+1
	  		print(f,c[n],s,"_",x,"_ALT",delim)
		else
		    n=n+1
	  		print(f,c[n],s,"_",x,delim)
		end	
	end
end
print(f,"\n" )

 

count = 0
v = "dmy"
while !eof(vc)
    v = next1(vc)
    count = count + 1
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
