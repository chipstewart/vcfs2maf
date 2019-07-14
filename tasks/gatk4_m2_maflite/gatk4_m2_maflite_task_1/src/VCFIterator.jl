export VCFIterator, getVersion, getFilename, getSampleNames, getINFOProperties, getFORMATProperties, getContigs, getFilters, getReference, next1, eof, close, reset
import Base.eof, Base.close, Base.reset
import GZip

#Main class for VCFiterator
#type VCFIterator
mutable struct VCFIterator
    version::AbstractString
    filename::AbstractString
    samples::Array{AbstractString, 1}
    infoTypes::Dict{AbstractString, Any}
    infoFlags::Array{AbstractString, 1}
    formatTypes::Dict{AbstractString, Any}
    contigs::Dict{AbstractString, Int64}
    filters::Dict{AbstractString, AbstractString}
    reference::AbstractString
    file 
end

#Stores INFO and FORMAT section properties
#type InfoField
mutable struct InfoField
    Number::AbstractString 
    Type::AbstractString
    Description::AbstractString
end

#Standard get methods
getVersion(vc::VCFIterator) = vc.version
getFilename(vc::VCFIterator) = vc.filename
getSampleNames(vc::VCFIterator) = vc.samples
getINFOProperties(vc::VCFIterator) = vc.infoTypes
getFORMATProperties(vc::VCFIterator) = vc.formatTypes
getContigs(vc::VCFIterator) = vc.contigs
getFilters(vc::VCFIterator) = vc.filters
getReference(vc::VCFIterator) = vc.reference

#Gets next variant.
function next1(vc::VCFIterator)
    if !eof(vc.file)
        line = readline(vc.file)
        if length(line)<2
            return false
        end
        return generateVariant(vc, line)
    else
        return false
    end
end

#Checks if it is eof.
function eof(vc::VCFIterator)
    eof(vc.file)
end

#Closes the parser
function close(vc::VCFIterator)
    close(vc.file)
end

#resets the parser to the first variant of the VCF file.
function reset(vc::VCFIterator)
    close(vc.file)
    vc.file = open(vc.filename)
    f2=splitext(filename)
    if (uppercase(f2[2])==".GZ")
        close(vc.file)
        f=GZip.open(vc.filename)
    end

    while !eof(vc.file)
        line = chomp(readline(vc.file))
        if line[1:2] != "##"
            break
        end
    end
end

#Parses info/format line
function readInfoLine!(typesdict::Dict{AbstractString, Any}, line::AbstractString, flaglist)
    #parse the line
    entries = split(line[2:end-1], ',')
    Id = split(entries[1], '=')[2]
    Number = split(entries[2], '=')[2]
    Type = split(entries[3], '=')[2]
    Description = split(entries[4], '=')[2]
    temp = InfoField(Number, Type, Description[2:end-1])
    
    #add to dict
    typesdict[Id] = temp
    
    if Type == "Flag"
        push!(flaglist, Id)
    end
end

#Parses contig line
function readContigLine!(contigs::Dict{AbstractString, Int64}, line::AbstractString)
    entries = split(line[2:end-1], ',')
    Id = split(entries[1], '=')[2]
    #println(entries[2])
    Length = parse(Int32,split(entries[2], '=')[2])
    contigs[Id] = Length
end

#Parses filter line
function readFilterLine!(filters::Dict{AbstractString, AbstractString}, line::AbstractString)
    entries = split(line[2:end-1], ',')
    Id = split(entries[1], '=')[2]
    Description = split(entries[2], '=')[2]
    filters[Id] = Description[2:end-1]
end

#Constructor for VCFIterator
function VCFIterator(filename::AbstractString)
    f=open(filename)
    f2=splitext(filename)
    if (uppercase(f2[2])==".GZ")
        close(f)
        f=GZip.open(filename)
    end
    vc = VCFIterator("VCFv4.1", filename, [], Dict{AbstractString, Any}(), [], Dict{AbstractString, Any}(), Dict{AbstractString, Int64}(), Dict{AbstractString, AbstractString}(), "", f)
    parseHeader(vc)
    
    #return at the end 
    vc 
end

#Parses through header lines
function parseHeader(vc::VCFIterator)
    while !eof(vc.file)
        line = chomp(readline(vc.file))
        if line[1:2] != "##"
            samples = split(line, "\t")[10:end]
            vc.samples = samples
            break
        end
        
        #read the info line
        if length(line) >= 6 && line[1:6]=="##INFO"
            line = line[8:end]
            readInfoLine!(vc.infoTypes, line, vc.infoFlags)
        elseif length(line) >= 8 && line[1:8]=="##FORMAT"
            line = line[10:end]
            readInfoLine!(vc.formatTypes, line, AbstractString[])
        elseif length(line) >= 8 && line[1:8]=="##contig"
            line = line[10:end]
            readContigLine!(vc.contigs, line)
        elseif length(line) >= 8 && line[1:8]=="##FILTER"
            line = line[10:end]
            readFilterLine!(vc.filters, line)
        elseif length(line) >= 11 && line[1:11]=="##reference"
            line = line[13:end]
            vc.reference = line
        end    
            
    end 
end

#Helper function that convertss x to float except when x is a period.
function floatField(x)
    if x == "."
        return "."
    else
        return parse(Float64,x)
    end
end

#Helper function that convertss x to int except when x is a period.
function intField(x)
    if x == "."
        return "."
    else
        #println(x)
        return parse(Int32,x)
    end
end

#Constructor for variant given a line of vcf file
function generateVariant(vc::VCFIterator, line::AbstractString)
    #remove newline and split by tabs. 
    line = split(chomp(line), "\t")
    
    #Assigining variables
    chr = convert(AbstractString, line[1])
    pos = intField(line[2])
    id_ = split(line[3])
    ref = convert(AbstractString, line[4])
    alt = split(line[5],",")
    if line[6] == "."
        qual = "."
    else
        qual = float(line[6])
    end
    filter = split(line[7],";")
    
    #creating INFO field
    information = Dict{AbstractString, Any}()
    tempINFO = split(line[8],";")
    
    #instantiate all flags to be false
    if tempINFO[1] != "."
        for i in 1:length(vc.infoFlags)
            information[vc.infoFlags[i]] = false
        end
        for i in 1:length(tempINFO)
            entry = tempINFO[i]
            id = split(entry,'=')[1]
            if ! haskey(vc.infoTypes, id) 
                continue
            end
            #check if it is flag or not
            contentType = (vc.infoTypes[id]).Type

            if contentType != "Flag"
                contents = split(entry,'=')[2]
                #parse through contents and convert types
                contents = split(contents, ",")

                if contentType == "Integer"
                    contents = map(intField, contents)
                elseif contentType == "Float"
                    contents = map(floatField, contents)
                elseif contentType != "String" && contentType != "Char"
                    println("INVALID TYPE",":\t",entry)
                end

                information[id] = contents
            else
                #set variable to true if the flag exists
                information[id] = true 
            end
        end
    end
    
    #creating FORMAT fields if they exist
    if length(line)>8
        format = Dict{AbstractString, Dict{AbstractString, Any}}()
        formatEntries = split(line[9],":")
        for i in 10:length(line)
            tempFORMAT = split(line[i],":")
            sample = vc.samples[i-9]
            tempFORMATdict = Dict{AbstractString, Any}()
            
            for j in 1:length(tempFORMAT)
                id = formatEntries[j]
                contents = split(tempFORMAT[j], ",")
                if ! haskey(vc.formatTypes, id) 
                    continue
                end
        
                contentType = vc.formatTypes[id].Type

                if contentType == "Integer"
                    contents = map(intField, contents)
                elseif contentType == "Float"
                    contents = map(floatField, contents)
                elseif contentType != "String" && contentType != "Char"
	            println("INVALID TYPE",":\t",entry)
                end
            
                tempFORMATdict[id] = contents
            end
            format[sample] = tempFORMATdict
        end
    else
        format = Dict{AbstractString,Dict{AbstractString,Any}}()
    end

    #construct a variant
    Variant(chr, pos, id_, ref, alt, qual, filter, information, format)
end
