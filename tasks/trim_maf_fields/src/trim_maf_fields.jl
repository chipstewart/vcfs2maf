#!/usr/local/bin/julia
#ARGS
# ARGS=["THCA-EM-A2CP-TP-NB", "THCA-EM-A2CP-TP-NB.annotated.maf",  "fields_to_remove.txt", "THCA-EM-A2CP-TP-NB.annotated.trim.maf"]

emptymaffields=split("Hugo_Symbol Entrez_Gene_Id  Center  NCBI_Build  Chromosome  Start_position  End_position    Strand  Variant_Classification  Variant_Type    Reference_Allele    Tumor_Seq_Allele1   Tumor_Seq_Allele2   dbSNP_RS    dbSNP_Val_Status    Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1    Tumor_Validation_Allele2    Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status Validation_Status   Mutation_Status Sequencing_Phase    Sequence_Source Validation_Method   Score   BAM_file    Sequencer   Tumor_Sample_UUID   Matched_Norm_Sample_UUID    Genome_Change   Annotation_Transcript   Transcript_Strand   Transcript_Exon Transcript_Position cDNA_Change Codon_Change    Protein_Change  ref_context gc_content  COSMIC_n_overlapping_mutations  M2  M2_DP   M2_ECNT M2_IN_PON   M2_NLOD M2_NLOD_full    M2_NORMAL_AD_ALT    M2_NORMAL_AD_ALT_full   M2_NORMAL_AD_REF    M2_NORMAL_AD_REF_full   M2_N_ART_LOD    M2_POP_AF   M2_P_GERMLINE   M2_RPA  M2_RU   M2_STR  M2_TLOD M2_TLOD_full    M2_TUMOR_AD_ALT M2_TUMOR_AD_ALT_full    M2_TUMOR_AD_REF M2_TUMOR_AD_REF_full    M2_TUMOR_AF M2_TUMOR_AF_full    M2_TUMOR_ALT_F1R2   M2_TUMOR_ALT_F2R1   M2_TUMOR_FOXOG  M2_TUMOR_FOXOG_full M2_TUMOR_MBQ    M2_TUMOR_MCL    M2_TUMOR_MFRL   M2_TUMOR_MMQ    M2_TUMOR_MPOS   M2_TUMOR_OBAMRC M2_TUMOR_REF_F1R2   M2_TUMOR_REF_F1R2_full  M2_TUMOR_REF_F2R1   M2_TUMOR_REF_F2R1_full  M2_full M2_n_alt_count  M2_n_ref_count  M2_normal_barcode   M2_t_alt_count  M2_t_lod_fstar  M2_t_lod_fstar_full M2_t_ref_count  M2_tumor_barcode    M2_tumor_f  M2_variant_type NALG    Production  Production_fstar_tumor_lod  Production_t_alt_count  Production_t_ref_count  Production_variant_type ccds_id gencode_transcript_name gencode_transcript_status   gencode_transcript_tags gencode_transcript_type gene_type   n_alt_count n_alt_count_full    n_ref_count n_ref_count_full    strelka strelka_NORMAL_ACGT_TIR_TOR strelka_NORMAL_DP   strelka_QS  strelka_TUMOR_ACGT_TIR_TOR  strelka_TUMOR_DP    strelka_n_alt_count strelka_n_ref_count strelka_normal_barcode  strelka_t_alt_count strelka_t_ref_count strelka_tumor_barcode   svaba   svaba_DBSNP svaba_LOD   svaba_MAPQ  svaba_NM    svaba_REPSEQ    svaba_SPAN  svaba_n_alt_count   svaba_n_ref_count   svaba_normal_CR svaba_normal_DP svaba_normal_LO svaba_normal_SR svaba_normal_barcode    svaba_t_alt_count   svaba_t_ref_count   svaba_tumor_CR  svaba_tumor_DP  svaba_tumor_LO  svaba_tumor_SR  svaba_tumor_barcode t_alt_count t_alt_count_full    t_ref_count t_ref_count_full")


using DataFrames

id=ARGS[1]
in_file=ARGS[2]
remove_file=ARGS[3]
out_file=ARGS[4]
#comment_char=ARGS[5]
#separator=ARGS[6]

f = open(remove_file);
remove = readlines(f)
close(f)

f = open(in_file);
skip = -1
x="#"
while x[1]=='#'
     x = readline(f)
     skip += 1
end
close(f)

#Strelka SNV vcf
df = readtable(in_file,separator='\t',skipstart = skip, nastrings=[""])
f = names(df)
r=map(x -> Symbol(chomp(x)),remove)

k=findin(r, f)
r=r[k]

delete!(df, r)

n=length(df[:Hugo_Symbol])
if (n<1)
    open(out_file, "w") do f
        writedlm(f, reshape(emptymaffields, 1, length(emptymaffields)), '\t')
    end
    quit()
end

isString(x::Number)=false
isString(x::DataArrays.NAtype)=false
isString(x::AbstractString)=true

maf=df
for c in names(maf)
    #println(c)
    if ~isString(maf[1,c])
        maf[c] = map(x -> string(x),maf[c])
    end
    k=find(map(x -> isna(x),maf[c]))
    maf[k,c]=""
end


for c in names(maf)
    println(c)

    k=find(map(x->x=="NA",maf[c]))
    maf[k,c] = ""
    k=find(map(x->uppercase(x)=="__UNKNOWN__",maf[c]))
    maf[k,c] = ""
end

open(out_file, "w") do f
    writedlm(f, reshape(names(maf), 1, length(names(maf))), '\t')
    writedlm(f, convert(Array,maf), '\t')
end


