#1 python
'''
Created on 18 Feb 2013=7

@author: stewart
'''
import csv
import argparse
import sys
import os
from quicksect import IntervalNode
verbose=False

if not (sys.version_info[0] == 2  and sys.version_info[1] in [7]):
    raise "Must use Python 2.7.x"


def parseOptions():
    description = '''
    Given a tsv file with headers on the first line, impose a genomic interval list filter.
    '''

    epilog= '''

        Writes subset of tsv file contained interval list for each sample

        Required columns in the input file (case sensitive):
             chromsome_field (Chromosome)
             position_field (Start_position)
                          
            '''
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-s','--sample_id', metavar='sample_id', type=str, help ='sample id.')
    parser.add_argument('-l','--intervalListFile', metavar='intervalListFile', type=str, help ='intervalList file... chr start end')
    parser.add_argument('-i','--inputFile', metavar='inputFile', type=str, help ='inputFile tsv.')
    parser.add_argument('-v','--verbose', action='store_true', help ='verbose messages.')
    parser.add_argument('-c','--chromsome_field',  metavar='chromsome_field', type=str, help ='chromsome field name in tsv.',default='Chromosome')
    parser.add_argument('-p','--position_field',  metavar='position_field', type=str, help ='position field name in tsv.',default='Start_position')
    parser.add_argument('-o','--output_area', metavar='output_area', type=str, help ='output area.',default='.')
    parser.add_argument('-f','--filter_stub', metavar='filter_stub', type=str, help ='filter_stubstring.',default='')
    parser.add_argument('-b','--blacklist', action='store_true', help ='interval list is a blacklist.')
    parser.add_argument('-C','--COMMENT_CHAR', metavar='COMMENT_CHAR', type=str, help ='COMMENT_CHAR at start of line.',default='#')
    args = parser.parse_args()

    return args

def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end) for x in out ]


class intervalList():
    """
    load maf into container object
    store header lines identified with # in first column
    dict with key Chromosome:Start_position:tumor_Sample_Barcode:Matched_Norm_Sample_Barcode
    """
    def __init__(self, fileName):

        self.N = 0
        self.headline=[]
        self.tree={}

        if not os.path.exists(fileName):
            return


        fh = open(fileName,'rt')
        a=-1
        for lines in fh:
            line=lines.strip()
            if line[0] in ['#','@']:
                self.headline.append(line)
                continue

            IL1 = line.split('\t')

            self.N=self.N+1

            chrom         = IL1[0]
            pstart        = IL1[1]
            pend          = IL1[2]

            p1=int(pstart)
            p2=int(pend)
            #chrom=chrom.replace("MT","M")
            #self.IL.append([chrom,p1,p2])

            start, end = p1-1, p2+1
            if not (chrom==a):
                if not (a==-1):
                    self.tree[a]=tree1
                    if (verbose):
                        print('scan %s  line %d chrom %s' % (fileName,self.N, a))
                a=chrom
                tree1 = IntervalNode( start, end )
            else: # build an interval tree from the rest of the data
                tree1 = tree1.insert( start, end )

        self.tree[a]=tree1
        fh.close()
        print "Loaded IntervalList:\t" + fileName
        print "events:\t" + str(self.N)



if __name__ == '__main__':

    args = parseOptions()
    id = args.sample_id
    intervalListFile = args.intervalListFile
    inputFile = args.inputFile
    CHROMSOME = args.chromsome_field
    POSITION = args.position_field
    verbose = args.verbose
    output = args.output_area
    filter_stub = args.filter_stub
    blacklist = args.blacklist
    COMMENT_CHAR = args.COMMENT_CHAR

    # load the input files
    IntervalList = intervalList(intervalListFile)

    if output is None:
        output = "."

    if not os.path.exists(output):
        os.makedirs(output)

    (p,inputFileName) = os.path.split(inputFile)
    (q,ext) = os.path.splitext(inputFile)

    print "id:\t" + id + "\n"
    print "input file:\t" + inputFile + "\n"
    print "interval list file:\t" + intervalListFile + "\n"


    output_Filename = output + "/" + id + filter_stub + ext

    print "output file:\t" + output_Filename + "\n"


    outputFileFP = file(output_Filename, 'w')

    # loop over call stats files
    fh = open(inputFile,'rt')
    nin=0
    npass=0
    nreject=0
    nhead = 0
    kC=-1
    kP=-1
    for lines in fh:
        line=lines.strip()
        if line[0] == COMMENT_CHAR:
            nhead = nhead + 1
            outputFileFP.write(line)
            continue
        if (nhead>0) and (nin<1):
            outputFileFP.write( COMMENT_CHAR + ' input ' + inputFile +  ' Filtered by interval list ' + intervalListFile + '\n')
            
        # first line after comment header lines should be field names 
        cs = line.split('\t')
        if (CHROMSOME in cs) and (POSITION in cs):
            outputFileFP.write(line+"\n")


            nC = cs.count(CHROMSOME)
            nP = cs.count(POSITION)
            if (nC*nP)<1:
                print "input missing required fields"
                if (nC < 1):
                    print  CHROMSOME
                    print line
                if (nP < 1):
                    print  POSITION
                    print line
                print "missing fields are allowed if input has no events"
                continue

            kC = cs.index(CHROMSOME)
            kP = cs.index(POSITION)


            if (kC<0) or (kP<0) :
                print "call stats contains incomplete records"
                if (kC < 0):
                    print  CHROMSOME
                    print line                    
                if (kP < 0):
                    print  POSITION
                    print line                    
                raise NameError('Cinput missing required fields ')
            continue

        if kC<0 or kP<1:
            print "input missing required fields"
            if (kC < 0):
                print  CHROMSOME
                print line
            if (kP < 0):
                print  POSITION
                print line
            raise NameError('Cinput missing required fields ')

        chrom         = cs[kC]
        position      = cs[kP]
        nin=nin+1

        #chrom=chrom.replace("MT","M")

        p1=int(position)
        p2=p1
        overap=[]
        if chrom in IntervalList.tree.keys():
            overlap = find(p1, p2 , IntervalList.tree[chrom])
        
        OK = len(overlap)>0
        
        if blacklist:
            OK = not OK


        if OK:
            outputFileFP.write(line + "\n")
            npass=npass+1
            if (verbose):
                print("pass:\t%d (%s:%s)\t of total %d" % (npass,chrom,position,nin))

        else:
            nreject=nreject+1
            if (verbose):
                print("reject:\t%d (%s:%s)\t of total %d" % (npass,chrom,position,nin))


    fh.close()

    outputFileFP.close()

    print "\n"
    print "input  #calls:\t" + str(nin)
    print "pass   #calls:\t" + str(npass)
    print "reject #calls:\t" + str(nreject) + "\n"
