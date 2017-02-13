export Variant, getChrom, getPos, getId, getRef, getAlt, getQual, getFilter, getInfo, getFormat

#Stores data for each variant
type Variant
    CHROM::String
    POS::Int64
    ID::Array{String,1}
    REF::String
    ALT::Array{String,1}
    #QUAL::Array{Float64, 1}
    QUAL
    FILTER::Array{String,1} 
    INFO::Dict{String, Any}
    FORMAT::Dict{String, Dict{String, Any}}
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

