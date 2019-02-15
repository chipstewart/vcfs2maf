module jlVCF

export Variant, getChrom, getPos, getId, getRef, getAlt, getQual, getFilter, getInfo, getFormat
export VCFIterator, getVersion, getFilename, getSampleNames, getINFOProperties, getFORMATProperties, getContigs, getFilters, getReference, next1, eof, close, reset

include("Variant.jl")
include("VCFIterator.jl") 

end 
