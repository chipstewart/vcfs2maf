#!/usr/local/bin/julia
# ARGS=["/Users/stewart/GoogleDrive/SU2C/FC/DF0051-N.hybrid_selection_metrics","MEAN_TARGET_COVERAGE"]
in_file=ARGS[1]
FIELD=ARGS[2]
out_file=string(FIELD , ".txt")

print(ARGS)

if !isfile(in_file)
    println("missing input maf ", in_file, " .")
end

f=open(in_file)
data1=""
while !eof(f)
    line = chomp(readline(f))
    if isempty(line)
      continue
    end
    if line[1] != '#'
      fields = split(line, "\t")
      data = split(chomp(readline(f)),"\t")
      k=find([contains(FIELD,i) for i in fields])
      data1=data[k]
      break
    end
close(f)
end
                 

open(out_file, "w") do fo
    writedlm(fo,data1)
end