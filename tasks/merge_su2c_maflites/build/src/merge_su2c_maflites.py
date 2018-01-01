from argparse import ArgumentParser, RawDescriptionHelpFormatter
#import csv
#import os
#import re
import sys
import numpy as np
import pandas as pd
from pandas import DataFrame

def parseOptions():
    epilog = """Merges four algs into one maflite"""
    desc = "merge input mafs ."
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)

    parser.add_argument("-P","--production_maf", type=str, help="Input production_maf")
    parser.add_argument("-2","--M2_maflite", type=str, help="Input  M2 maflite")
    parser.add_argument("-S","--SvABA_maflite", type=str, help="Input SvABA_maflite")
    parser.add_argument("-s","--strelka_maflite", type=str, help="Input strelka maflite")
    parser.add_argument("-i","--id", type=str, help="id stub prepended to output files.  <id>.merged.maflite.tsv")
    parser.add_argument("-t","--tumor_id", type=str, help="tumor sample barcode")
    parser.add_argument("-n","--normal_id", type=str, help="normal sample barcode")
    args = parser.parse_args()
    return args

def fill_in_end_and_var_type(row):
    if ((np.isnan(row['ref_allele']==True)) or (np.isnan(row['alt_allele']==True)) or (row['ref_allele'] is None) or (row['alt_allele'] is None)):
        end_pos = float('nan')
        var='nan'
        #print(var, end_pos)
        return end_pos,var #pd.DataFrame({'a':end_pos, 'b':'err'})
    r = str(row['ref_allele']).replace('-','')
    a = str(row['alt_allele']).replace('-','')
    if len(r) < len(a):
        end_pos = row['start']+1
        var = 'INS'
        #print(var, end_pos)
        return end_pos, var #pd.DataFrame({'a':end_pos, 'b':'INS'})#[end_pos, 'INS']
    if len(r) == len(a):
        end_pos = row['start']+len(str(row['alt_allele']))-1
        var = 'SNP'
        #print(var, end_pos)
        return end_pos, var
    var = 'DEL'
    end_pos = row['start']+len(str(row['ref_allele']))-1
    #print(var, end_pos)
    return end_pos, var


def svaba_alt_update(row):
    if(row['svaba_variant_type']=='DEL'):
        # take out the first character
        # Chip: do i need to change the start as well
        alt_update = row['alt_allele'][1:]
        if (alt_update ==''):
            alt_update='-'
        ref_update = row['ref_allele'][1:]
        if (ref_update ==''):
            ref_update='-'
        return alt_update,ref_update
    else:
        alt_update = row['alt_allele']
        ref_update = row['ref_allele']
        return alt_update,ref_update


# delete row when the value is 1
def M2_fill_in_multi_allele_var(row):
    # print("Start")
    # update: M2_val_multi','M2_val_multi_index and alt_allele
    alt_alleles = row['alt_allele']
    if alt_alleles is None:
        val_multi = str('nan')
        val_multi_index = str('nan')
        # print("return all null")
        delete_row = 0
        return val_multi, val_multi_index, alt_alleles, delete_row

    list_alleles = alt_alleles.split(':')

    if (len(list_alleles) == 1):  # only one allele
        # print("list alleles is 1")
        val_multi = list_alleles
        val_multi_index = '0'
        delete_row = 0
        # print(val_multi,val_multi_index,alt_alleles)
        return val_multi, val_multi_index, alt_alleles, delete_row

    alt_tumor_af = row['M2_TUMOR_AF']
    list_tumor_af = alt_tumor_af.split(':')

    if (len(list_alleles) != len(list_tumor_af)):
        val_multi = str('nan')
        val_multi_index = str('nan')
        delete_row = 1
        # print(val_multi,val_multi_index,alt_alleles)
        return val_multi, val_multi_index, alt_alleles, delete_row

    # assuming that the allele_fraction is allighned with the allele column order, use the allele fraction to select the index
    # of the allele with the highest allele fraction
    val_multi = alt_alleles
    alt_tumor_af_float = list(map(float, list_tumor_af))
    index_of_max_af = alt_tumor_af_float.index(max(alt_tumor_af_float))
    allele = list_alleles[index_of_max_af]
    delete_row = 0
    # print(val_multi,index_of_max_af,allele)
    return val_multi, index_of_max_af, allele, delete_row

def M2(input_tsv):
    # "SU2CLC-MGH-1048_1.fixed.m2.pass.maflite.tsv"
    # key - ['chr','start','end','alt_allele','ref_allele']
    columns_maf_lite_M2 = ['chr',
                           'start',
                           'end',
                           'ref_allele',
                           'alt_allele',
                           'M2_variant_type',
                           'M2_val_multi',
                           'M2_val_multi_index',
                           'M2_DP',
                           'M2_POP_AF',
                           'M2_IN_PON',
                           'M2_RPA',
                           'M2_P_GERMLINE',
                           'M2_STR',
                           'M2_RU',
                           'M2_NLOD',
                           'M2_ECNT',
                           'M2_N_ART_LOD',
                           'M2_TLOD',
                           'M2_TUMOR_MMQ',
                           'M2_TUMOR_MBQ',
                           'M2_TUMOR_AF',
                           'M2_TUMOR_FOXOG',
                           'M2_TUMOR_MPOS',
                           'M2_TUMOR_ALT_F2R1',
                           'M2_TUMOR_AD_REF',
                           'M2_TUMOR_AD_ALT',
                           'M2_TUMOR_REF_F2R1',
                           'M2_TUMOR_MCL',
                           'M2_TUMOR_MFRL',
                           'M2_TUMOR_OBAMRC',
                           'M2_TUMOR_REF_F1R2',
                           'M2_TUMOR_ALT_F1R2',
                           'M2_NORMAL_AD_REF',
                           'M2_NORMAL_AD_ALT',
                           'M2_tumor_barcode',
                           'M2_normal_barcode',
                           'M2_n_alt_count',
                           'M2_n_ref_count',
                           'M2_t_lod_fstar',
                           'M2_t_alt_count',
                           'M2_t_ref_count',
                           'M2_multi_allele_variant',
                           'M2',
                           'delete_row']
    # delete_row is used as a marker and will be taken out before saving to CSV

    df_M2 = DataFrame.from_csv(input_tsv, sep="\t", index_col=None)
    df_M2 = df_M2[df_M2.TUMOR_AD_ALT != 0]
    df_maf_lite_M2 = pd.DataFrame(index=None, columns=columns_maf_lite_M2)

    # START: MUST HAVE IN MAFLITE
    df_maf_lite_M2['chr'] = df_M2['chr']
    # +1 - Chip, I have removed the +1 because there is no variant type
    df_maf_lite_M2['start'] = df_M2['start']
    df_maf_lite_M2['ref_allele'] = df_M2['ref_allele']
    df_maf_lite_M2['alt_allele'] = df_M2['alt_allele']
    df_maf_lite_M2['M2_n_alt_count'] = df_M2['n_alt_count']
    df_maf_lite_M2['M2_n_ref_count'] = df_M2['NORMAL_AD_REF']
    df_maf_lite_M2['M2_tumor_f'] = df_M2['TUMOR_AF']
    df_maf_lite_M2['M2_t_alt_count'] = df_M2['TUMOR_AD_ALT']
    df_maf_lite_M2['M2_t_ref_count'] = df_M2['TUMOR_AD_REF']
    # END: MUST HAVE IN MAFLITE

    df_maf_lite_M2['M2_DP'] = df_M2['DP']
    df_maf_lite_M2['M2_POP_AF'] = df_M2['POP_AF']
    df_maf_lite_M2['M2_IN_PON'] = df_M2['IN_PON']
    df_maf_lite_M2['M2_RPA'] = df_M2['RPA']
    df_maf_lite_M2['M2_P_GERMLINE'] = df_M2['P_GERMLINE']
    df_maf_lite_M2['M2_STR'] = df_M2['STR']
    df_maf_lite_M2['M2_RU'] = df_M2['RU']
    df_maf_lite_M2['M2_NLOD'] = df_M2['NLOD']
    df_maf_lite_M2['M2_ECNT'] = df_M2['ECNT']
    df_maf_lite_M2['M2_N_ART_LOD'] = df_M2['N_ART_LOD']
    # df_maf_lite['TLOD']=df['TLOD']
    df_maf_lite_M2['M2_TUMOR_MMQ'] = df_M2['TUMOR_MMQ']
    df_maf_lite_M2['M2_TUMOR_MBQ'] = df_M2['TUMOR_MBQ']
    df_maf_lite_M2['M2_TUMOR_AF'] = df_M2['TUMOR_AF']
    df_maf_lite_M2['M2_TUMOR_FOXOG'] = df_M2['TUMOR_FOXOG']
    df_maf_lite_M2['M2_TUMOR_MPOS'] = df_M2['TUMOR_MPOS']
    df_maf_lite_M2['M2_TUMOR_ALT_F2R1'] = df_M2['TUMOR_ALT_F2R1']
    df_maf_lite_M2['M2_TUMOR_AD_REF'] = df_M2['TUMOR_AD_REF']
    df_maf_lite_M2['M2_TUMOR_AD_ALT'] = df_M2['TUMOR_AD_ALT']
    df_maf_lite_M2['M2_TUMOR_REF_F2R1'] = df_M2['TUMOR_REF_F2R1']
    df_maf_lite_M2['M2_TUMOR_MCL'] = df_M2['TUMOR_MCL']
    df_maf_lite_M2['M2_TUMOR_MFRL'] = df_M2['TUMOR_MFRL']
    df_maf_lite_M2['M2_TUMOR_OBAMRC'] = df_M2['TUMOR_OBAMRC']
    df_maf_lite_M2['M2_TUMOR_REF_F1R2'] = df_M2['TUMOR_REF_F1R2']
    df_maf_lite_M2['M2_TUMOR_ALT_F1R2'] = df_M2['TUMOR_ALT_F1R2']
    df_maf_lite_M2['M2_NORMAL_AD_REF'] = df_M2['NORMAL_AD_REF']
    df_maf_lite_M2['M2_NORMAL_AD_ALT'] = df_M2['NORMAL_AD_ALT']
    df_maf_lite_M2['M2_t_lod_fstar'] = df_M2['TLOD']
    df_maf_lite_M2['M2_tumor_barcode'] = df_M2['tumor_barcode']
    df_maf_lite_M2['M2_normal_barcode'] = df_M2['normal_barcode']


    # Init end based on the start and the ref_Allele
    df_maf_lite_M2[['end', 'M2_variant_type']] = pd.Series.tolist(
        df_maf_lite_M2.apply(fill_in_end_and_var_type, axis=1, reduce=True))
    # select the allele with most evidence
    df_maf_lite_M2[['M2_val_multi', 'M2_val_multi_index', 'alt_allele', 'delete_row']] = pd.Series.tolist(
        df_maf_lite_M2.apply(M2_fill_in_multi_allele_var, axis=1))

    df_maf_lite_M2[df_maf_lite_M2['delete_row'] == 1]


    # delete rows that had delete_row =1, means that the number of alleles is different than their coverage
    df_maf_lite_M2 = df_maf_lite_M2[df_maf_lite_M2['delete_row'] == 0]
    # df_maf_lite_M2[df_maf_lite_M2.M2_val_multi_index != '1']

    # take out the column delete_row
    df_maf_lite_M2 = df_maf_lite_M2.drop('delete_row', 1)
    df_maf_lite_M2 = df_maf_lite_M2.drop('M2_val_multi', 1)
    df_maf_lite_M2 = df_maf_lite_M2.drop('M2_val_multi_index', 1)
    df_maf_lite_M2 = df_maf_lite_M2.drop('M2_multi_allele_variant', 1)

    # init M2 column with 1
    df_maf_lite_M2['M2'] = 1

    return df_maf_lite_M2

def SvABA(input_tsv):
    # "SU2CLC-MGH-1048_1.SvABA.INDEL.tsv"
    df_maf_lite_svaba = DataFrame.from_csv(input_tsv, sep="\t", index_col=None)
    df_maf_lite_svaba.rename(columns={'REPSEQ': 'svaba_REPSEQ',
                                      'NM': 'svaba_NM',
                                      'MAPQ': 'svaba_MAPQ',
                                      'DBSNP': 'svaba_DBSNP',
                                      'LOD': 'svaba_LOD',
                                      'SPAN': 'svaba_SPAN',
                                      'normal_CR': 'svaba_normal_CR',
                                      'normal_DP': 'svaba_normal_DP',
                                      'normal_LO': 'svaba_normal_LO',
                                      'n_alt_count': 'svaba_n_alt_count',
                                      'normal_SR': 'svaba_normal_SR',
                                      'tumor_CR': 'svaba_tumor_CR',
                                      'tumor_DP': 'svaba_tumor_DP',
                                      'tumor_LO': 'svaba_tumor_LO',
                                      't_alt_count': 'svaba_t_alt_count',
                                      'tumor_SR': 'svaba_tumor_SR',
                                      'n_ref_count': 'svaba_n_ref_count',
                                      't_ref_count': 'svaba_t_ref_count',
                                      'tumor_barcode': 'svaba_tumor_barcode',
                                      'normal_barcode': 'svaba_normal_barcode',
                                        }, inplace=True)

    df_maf_lite_svaba['svaba'] = 1
    df_maf_lite_svaba = df_maf_lite_svaba.drop('judgement', 1)
    df_maf_lite_svaba = df_maf_lite_svaba.drop('build', 1)
    # save to Svaba maflite file
    # df_maf_lite_svaba.to_csv("SVABA_MAF_LITE_SU2CLC-MGH-1048_1.maf", sep="\t", index=None)

    return df_maf_lite_svaba

def SvABA0(input_tsv):
    # "SU2CLC-MGH-1048_1.SvABA.INDEL.tsv"
    df_svaba = DataFrame.from_csv(input_tsv, sep="\t", index_col=None)

    # create a MAFLITE for Svaba
    # Chip: do we need a different headers for Svaba? or the same like M2
    # TODO: add these names to the table: 'svaba_t_alt_count', 'svaba_t_ref_count', 'svaba_n_alt_count','svaba_n_ref_count'
    columns_maf_lite_svaba = ['chr', 'start', 'end', 'svaba_variant_type', 'ref_allele', 'alt_allele', 'svaba_IN_PON',
                              'svaba_t_alt_count', 'svaba_t_ref_count', 'svaba_n_alt_count', 'svaba_n_ref_count',
                              'svaba_multi_allele_variant', 'svaba']
    df_maf_lite_svaba = pd.DataFrame(index=None, columns=columns_maf_lite_svaba)

    # ['CHRO', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'REPSEQ', 'NM',
    #       'BX', 'PON', 'READNAMES', 'SCTG', 'MAPQ', 'DBSNP', 'LOD', 'SPAN',
    # START: MUST HAVE IN MAFLITE
    import pandas as pd
    df_maf_lite_svaba['chr'] = df_svaba['CHRO']
    df_maf_lite_svaba['start'] = df_svaba['POS'] + 1
    df_maf_lite_svaba['ref_allele'] = df_svaba['REF']
    df_maf_lite_svaba['alt_allele'] = df_svaba['ALT']
    df_maf_lite_svaba['svaba_IN_PON'] = df_svaba['PON']

    # update the string to be the next one
    df_maf_lite_svaba[['end', 'svaba_variant_type']] = pd.Series.tolist(
        df_maf_lite_svaba.apply(fill_in_end_and_var_type, axis=1))
    # df_maf_lite_svaba[['end','svaba_variant_type']] = pd.Series.to_frame(df_maf_lite_svaba.apply(fill_in_end_and_var_type, axis=1))
    # df_maf_lite_svaba.head()

    # relies on del
    df_maf_lite_svaba[['alt_allele', 'ref_allele']] = pd.Series.tolist(
        df_maf_lite_svaba.apply(svaba_alt_update, axis=1))
    df_maf_lite_svaba.head()

    df_maf_lite_svaba['alt_allele'][0] == ''

    l = df_svaba.columns
    l = list(l)
    svaba_bam_AD_T = None
    svaba_bam_AD_N = None
    svaba_bam_DP_T = None
    svaba_bam_DP_N = None

    svaba_bam_AD = [s for s in l if "bam_AD" in s]
    svaba_bam_DP = [s for s in l if "bam_DP" in s]

    print("len(svaba_bam_AD) ")
    print(len(svaba_bam_AD))
    print(svaba_bam_AD)

    print("len(svaba_bam_DP) ")
    print(len(svaba_bam_DP))
    print(svaba_bam_DP)

    if (len(svaba_bam_AD) != 2):
        print("error bam_AD")

    if (len(svaba_bam_DP) != 2):
        print("error bam_DP")

    if "TM" in svaba_bam_DP[0]:
        print("if")
        # df.rename(columns={'Leader': 'Commander'}, inplace=True)
        df_svaba.rename(columns={svaba_bam_DP[0]: 'svaba_bam_DP_T', svaba_bam_DP[1]: 'svaba_bam_DP_N'}, inplace=True)
    else:
        print("else")
        df_svaba.rename(columns={svaba_bam_DP[0]: 'svaba_bam_DP_N', svaba_bam_DP[1]: 'svaba_bam_DP_T'}, inplace=True)

    if "TM" in svaba_bam_AD[0]:
        df_svaba.rename(columns={svaba_bam_AD[0]: 'svaba_bam_AD_T', svaba_bam_AD[1]: 'svaba_bam_AD_N'}, inplace=True)
    else:
        df_svaba.rename(columns={svaba_bam_AD[0]: 'svaba_bam_AD_N', svaba_bam_AD[1]: 'svaba_bam_AD_T'}, inplace=True)

    # TODO: check that update all sample rows with the right values
    df_maf_lite_svaba['svaba_t_alt_count'] = df_svaba['svaba_bam_AD_T']
    df_maf_lite_svaba['svaba_t_ref_count'] = df_svaba['svaba_bam_DP_T'] - df_svaba['svaba_bam_AD_T']
    df_maf_lite_svaba['svaba_n_alt_count'] = df_svaba['svaba_bam_AD_N']
    df_maf_lite_svaba['svaba_n_ref_count'] = df_svaba['svaba_bam_DP_N'] - df_svaba['svaba_bam_AD_N']
    df_maf_lite_svaba['svaba'] = 1

    # save to Svaba maflite file
    #df_maf_lite_svaba.to_csv("SVABA_MAF_LITE_SU2CLC-MGH-1048_1.maf", sep="\t", index=None)

    return df_maf_lite_svaba

def strelka(input_tsv):
    # df_maf_lite_strelka = DataFrame.from_csv("SU2CLC-MGH-1048_1.Strelka_maflite.tsv", sep="\t",index_col=None)
    df_maf_lite_strelka = DataFrame.from_csv(input_tsv, sep="\t", index_col=None)
    # add rename
    df_maf_lite_strelka.rename(columns={'n_ref_count': 'strelka_n_ref_count',
                                        'NORMAL_DP': 'strelka_NORMAL_DP',
                                        'TUMOR_DP': 'strelka_TUMOR_DP',
                                        'n_alt_count': 'strelka_n_alt_count',
                                        't_ref_count': 'strelka_t_ref_count',
                                        't_alt_count': 'strelka_t_alt_count',
                                        'QS': 'strelka_QS',
                                        'NORMAL_ACGT_TIR_TOR': 'strelka_NORMAL_ACGT_TIR_TOR',
                                        'TUMOR_ACGT_TIR_TOR': 'strelka_TUMOR_ACGT_TIR_TOR',
                                        'tumor_barcode': 'strelka_tumor_barcode',
                                        'normal_barcode': 'strelka_normal_barcode',
                                        }, inplace=True)

    df_maf_lite_strelka['strelka'] = 1
    df_maf_lite_strelka = df_maf_lite_strelka.drop('judgement', 1)
    df_maf_lite_strelka = df_maf_lite_strelka.drop('build', 1)
    return df_maf_lite_strelka

def Production(input_maf):
    df_maf_production = DataFrame.from_csv(input_maf, sep="\t", index_col=None)
    df_maf_production.head(2)

    # columns_maf_lite_svaba = ['chr','start','end','svaba_variant_type','ref_allele','alt_allele','IN_PON']
    columns_maf_lite_production = ['chr', 'start', 'ref_allele', 'alt_allele', 'end', 'Production_t_alt_count',
                                   'Production_t_ref_count', 'Production_variant_type', 'Production_fstar_tumor_lod',
                                   'Production']
    # create a production maflite
    df_maf_lite_production = pd.DataFrame(index=None, columns=columns_maf_lite_production)
    # START: MUST HAVE IN MAFLITE
    df_maf_lite_production['chr'] = df_maf_production['Chromosome']
    df_maf_lite_production['start'] = df_maf_production['Start_position']
    df_maf_lite_production['ref_allele'] = df_maf_production['Reference_Allele']
    df_maf_lite_production['alt_allele'] = df_maf_production['Tumor_Seq_Allele2']
    df_maf_lite_production['end'] = df_maf_production['End_position']
    #df_maf_lite_production['prod_tumor_f'] = df_maf_production['tumor_f']
    df_maf_lite_production['Production_t_alt_count'] = df_maf_production['t_alt_count']
    df_maf_lite_production['Production_t_ref_count'] = df_maf_production['t_ref_count']
    # END: MUST HAVE IN MAFLITE
    df_maf_lite_production['Production_fstar_tumor_lod'] = df_maf_production['fstar_tumor_lod']
    df_maf_lite_production['Production'] = 1
    df_maf_lite_production.head(2)
    return df_maf_lite_production

def merge(df_maf_lite_M2, df_maf_lite_production,df_maf_lite_svaba,df_maf_lite_strelka):

    # TODO: check if is being sorted and how we can work on the indel - range
    # TODO: check how to take out the suffix
    M2_Prod_Merge = pd.merge(df_maf_lite_M2, df_maf_lite_production,
                             on=['chr', 'start', 'end', 'alt_allele', 'ref_allele'], how='outer', sort=True,
                             suffixes=('_M2', '_prod'), indicator=False)
    print(M2_Prod_Merge.head(2))
    Svaba_Stelka_Merge = pd.merge(df_maf_lite_svaba, df_maf_lite_strelka,
                                  on=['chr', 'start', 'end', 'alt_allele', 'ref_allele'], how='outer', sort=True,
                                  suffixes=('_svaba', '_strelka'), indicator=False)
    print(Svaba_Stelka_Merge.head(2))
    print("*************************")
    # Maybe i need to join
    all_merged = pd.merge(M2_Prod_Merge, Svaba_Stelka_Merge, on=['chr', 'start', 'end', 'alt_allele', 'ref_allele'],
                          how='outer', sort=True, suffixes=('_M2_Prod_Merge', '_Svaba_Stelka_Merge'), indicator=False)
    # print(all_merged.head(2))
    # now i can replace the left_only with M2 and the right only with M2, rename the table and generate another merge for the other two and then merge betweent he results
    # after that , check for 'start', 'end' & variant type (cechk that i have the variant type on all)
    #
    all_merged.loc[all_merged['Production'].isnull(), 'Production'] = 0
    all_merged.loc[all_merged['M2'].isnull(), 'M2'] = 0
    all_merged.loc[all_merged['svaba'].isnull(), 'svaba'] = 0
    all_merged.loc[all_merged['strelka'].isnull(), 'strelka'] = 0

    all_merged['NALG'] = all_merged['Production']+all_merged['M2']+all_merged['svaba']+all_merged['strelka']
    #
    #  allele counts from production first (only t counts, n counts from strelka)
    all_merged['t_alt_count'] =  all_merged['Production_t_alt_count']
    all_merged['t_ref_count'] =  all_merged['Production_t_ref_count']
    all_merged['n_alt_count'] =  all_merged['strelka_n_alt_count']
    all_merged['n_ref_count'] =  all_merged['strelka_n_ref_count']

    # then Strelka
    if (all_merged['t_alt_count'].isnull().values.any()):
        k = all_merged.index[all_merged['t_alt_count'].isnull()]
        for i in k:
            all_merged.loc[i, 't_alt_count'] = all_merged.loc[i, 'strelka_t_alt_count']
        k = all_merged.index[all_merged['t_ref_count'].isnull()]
        for i in k:
            all_merged.loc[i, 't_ref_count'] = all_merged.loc[i, 'strelka_t_ref_count']
       #M2
        k = all_merged.index[all_merged['t_alt_count'].isnull()]
        for i in k:
            all_merged.loc[i, 't_alt_count'] = all_merged.loc[i, 'M2_t_alt_count']
        k = all_merged.index[all_merged['t_ref_count'].isnull()]
        for i in k:
            all_merged.loc[i, 't_ref_count'] = all_merged.loc[i, 'M2_t_ref_count']
        k = all_merged.index[all_merged['n_alt_count'].isnull()]
        for i in k:
            all_merged.loc[i, 'n_alt_count'] = all_merged.loc[i, 'M2_n_alt_count']
        k = all_merged.index[all_merged['n_ref_count'].isnull()]
        for i in k:
            all_merged.loc[i, 'n_ref_count'] = all_merged.loc[i, 'M2_n_ref_count']

        #svaba
        k = all_merged.index[all_merged['t_alt_count'].isnull()]
        for i in k:
            all_merged.loc[i, 't_alt_count'] = all_merged.loc[i, 'svaba_t_alt_count']
        k = all_merged.index[all_merged['t_ref_count'].isnull()]
        for i in k:
            all_merged.loc[i, 't_ref_count'] = all_merged.loc[i, 'svaba_t_ref_count']
        k = all_merged.index[all_merged['n_alt_count'].isnull()]
        for i in k:
            all_merged.loc[i, 'n_alt_count'] = all_merged.loc[i, 'svaba_n_alt_count']
        k = all_merged.index[all_merged['n_ref_count'].isnull()]
        for i in k:
            all_merged.loc[i, 'n_ref_count'] = all_merged.loc[i, 'svaba_n_ref_count']

    if (all_merged['t_alt_count'].isnull().values.any()):
        print('missing alt counts')

    return all_merged


def main():
    args = parseOptions()

    production_maf = args.production_maf
    M2_maflite = args.M2_maflite
    SvABA_maflite = args.SvABA_maflite
    strelka_maflite = args.strelka_maflite
    id = args.id
    tumor_id = args.tumor_id
    normal_id = args.normal_id

    df_maf_lite_M2 = M2(M2_maflite)
    df_maf_lite_M2.to_csv(id+".M2.maflite.tsv", sep="\t", index=None)

    df_maf_lite_svaba = SvABA(SvABA_maflite)
    df_maf_lite_svaba.to_csv(id+".svaba.maflite.tsv", sep="\t", index=None)

    df_maf_lite_strelka = strelka(strelka_maflite)
    df_maf_lite_strelka.to_csv(id+".strelka.maflite.tsv", sep="\t", index=None)

    df_maf_lite_production = Production(production_maf)
    df_maf_lite_production.to_csv(id+".production.maflite.tsv", sep="\t", index=None)

    df_maf_lite_merge=merge(df_maf_lite_M2, df_maf_lite_production,df_maf_lite_svaba,df_maf_lite_strelka)
    df_maf_lite_merge['tumor_barcode']=tumor_id
    df_maf_lite_merge['normal_barcode']=normal_id
    df_maf_lite_merge.to_csv(id+".merge.maflite.tsv", sep="\t", index=None)

    return 0

if __name__ == "__main__":
    sys.exit(main())

