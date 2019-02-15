# module Variant

export Variant, getChrom, getPos, getId, getRef, getAlt, getQual, getFilter, getInfo, getFormat

#Stores data for each variant
#type Variant
mutable struct  Variant
    CHROM::AbstractString
    POS::Int64
    ID::Array{AbstractString,1}
    REF::AbstractString
    ALT::Array{AbstractString,1}
    #QUAL::Array{Float64, 1}
    QUAL
    FILTER::Array{AbstractString,1} 
    INFO::Dict{AbstractString, Any}
    FORMAT::Dict{AbstractString, Dict{AbstractString, Any}}
end

#Get methods for variants
getChrom(v::Variant) = v.CHROM
getPos(v::Variant) = v.POS
getId(v::Variant) = v.ID
getRef(v::Variant) = v.REF
getAlt(v::Variant) = v.ALT
getQual(v::Variant) = v.QUAL
getFilter(v::Variant) = v.FILTER
getInfo(v::Variant) = v.INFO
getFormat(v::Variant) = v.FORMAT

#end

