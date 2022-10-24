#! /usr/bin/env python3
# coding: utf-8

import sys, getopt
import csv
import re
import logging
import os
from Bio import SeqIO
from Bio.Seq import Seq


def revcomp(seq):
    bioseq = Seq(seq)
    return str(bioseq.reverse_complement())

def main(argv):
    minyear=0
    maxyear=10000
    try:
        opts, args = getopt.getopt(argv,"hi:r:")
    except getopt.GetoptError:
        print('countMutations.py -i <alignment file> -r <reference file>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('countMutations.py -i <alignment file> -r <reference file>')
            sys.exit(0)
        elif opt in ("-i"):
            alignfile = arg
        elif opt in ("-r"):
            reffile = arg
            
    logging.info(f'Alignment file is {alignfile}')
    logging.info(f'Reference File id {reffile}')

    refseq=""
    # Parse reference sequence
    with open(reffile) as file:
        for record in SeqIO.parse(file, "fasta"):
            refseq=record.seq
    
    print(f'Sample\tAnnot1\tAnnot\tType\tPosition\tRef\tComp\tTrinucSeq1-TrinucSeq2\tStrand\tRefContext')
    with open(alignfile) as file:
        for record in SeqIO.parse(file, "fasta"):
            compname=record.description
            compseq=record.seq
            annotFull=compname.split("-")[-1]
            annot=re.sub('\d','',annotFull)
            
            if compseq == "" :
                logging.error(f'Sequence {compname} is null')
                sys.exit(1)
            
            if len(compseq) != len(refseq) :
                logging.error(f'Comp and reference sequences have different lengths')
                sys.exit(1)

            
                
            for c in range(0, len(compseq)):
                
                if refseq[c] != compseq[c] and compseq[c] != 'N':
                    type="Others"
                    strand="+"
                    refcontext=str(refseq[c-1:c+2])
                    if refseq[c]=='C' and compseq[c]=="T":
                        type="C-T"
                    elif refseq[c]=="G" and compseq[c]=="A":
                        strand='-'
                        refcontext=revcomp(refcontext)
                        type="G-A"

                    print(f'{record.description}\t{annotFull}\t{annot}\t{type}\t{c}\t{refseq[c]}\t{compseq[c]}\t{refseq[c-1:c+2]}-{compseq[c-1:c+2]}\t{strand}\t{refcontext}')
        
if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
