#!/usr/local/bin/julia
# ARGS=["/Users/stewart/GoogleDrive/SU2C/FC/DF0051-N.hybrid_selection_metrics","MEAN_TARGET_COVERAGE"]
in_file=ARGS[1]
FIELD=ARGS[2]
out_file="out.txt"

print(ARGS)

if !isfile(in_file)
    println("missing input file ", in_file, " .")
end

f=open(in_file)
data1=""
line=""
data=[]
while !eof(f)
    line = chomp(readline(f))
    println(line)
    if isempty(line)
      continue
    end
    if line[1] != '#'
      fields = split(line, "\t")
      data = split(chomp(readline(f)),"\t")
      k=find([contains(FIELD,i) for i in fields])
      data1=data[k]
      println(data1)
      break
    end
end
close(f)
          
println(line)
println(data)
println(data1)


open(out_file, "w") do fo
    writedlm(fo,data1)
end