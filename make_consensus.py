#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:47:22 2020

@author: avierstr
make the consensus sequence of sequences in a file

For c implementation of Levenshtein:
    https://pypi.org/project/python-Levenshtein/
    sudo apt-get install python3-dev
    (needed: python3-setuptools, python3-pip, python3-wheel)
    python3 -m pip install python-Levenshtein
"""

from Bio import SeqIO
from Levenshtein import distance as l_distance
from Levenshtein import median as l_median
import argparse

#==============================================================================    
def distance(A1,A2):  # calculate the similarity of 2 sequences
    if len(A1)*1.05 < len(A2): # if A1 is much shorter than A2, it influences the similarity
        idenprev = 0
        templist = []
        for a in range(0, len(A2) - len(A1),2): 
            e = int(len(A1)*1.015) # leave room for inserts
            distance = l_distance(A1, A2[a:e])
            iden = round(1 - distance/len(A2[a:e]),3)
            if iden >= idenprev:
                templist.append([a,iden])
                idenprev = iden
            else:
                break
        templist.sort(key=lambda x: x[1]) # sort according to value
        a, iden = templist[-1] # get best value
    else:
        distance = l_distance(A1, A2)
        iden = round(1 - distance/len(A2),3)
        
    return iden
#==============================================================================
def get_arguments():
    
    def range_limited_float_type(arg):
        """ Type function for argparse - a float within some predefined bounds """
        try:
            f = float(arg)
        except ValueError:    
            raise argparse.ArgumentTypeError("Must be a floating point number")
        if f < 50 or f > 100:
            raise argparse.ArgumentTypeError("Argument must be > " + str(50.0) + " and < " + str(100.0))
        return f

    parser = argparse.ArgumentParser(description='make_consensus: makes the consensus of reads in a file.' )
    parser.add_argument('-i', '--input', required=True, 
                        help='Input file in fastq of fasta format')
    parser.add_argument('-s', '--similar', type = range_limited_float_type, required=False, default=85.0 ,
                        help='Similarity to add a read to the consensus (value between 50 and 100). Default=85.0')

    args = parser.parse_args()
    return args
#==============================================================================
def compl_reverse(self):
    inp  = 'ATCG' # translate table for complement
    outp = 'TAGC'
    complement = ''.maketrans(inp, outp)
    R = (self[::-1]).translate(complement)  # complement reverse
    return R
#==============================================================================
def consensus_direction(consensuslist): #  check if all sequences are in the same direction (F or R)
    similar = args.similar/100
    consensusset = set()
    for x in consensuslist[0:1]:
        for y in consensuslist[1:]:
            iden = distance(x,y)
            idenR = distance(x,compl_reverse(y))
            if iden > idenR and iden >= similar:
                consensusset.update([x, y])
            else:
                if idenR >= similar:
                    consensusset.update([x, compl_reverse(y)])
    consensuslist = list(consensusset)
    return consensuslist
#==============================================================================
def read_file(infile):
    similar = args.similar
    with open(infile, 'r') as inf: # check the fileformat
        line = inf.readline()
        if line[0] == '>':
            fileformat = 'fasta'
        elif line[0] == '@':
            fileformat = 'fastq'
            
    comparelist = []
    inputfile = open(infile, "r")
    for record in SeqIO.parse(inputfile, fileformat):
        comparelist.append(str(record.seq))
    inputfile.close()
    
    comparelist2 = consensus_direction(comparelist) # get all seq in same direction
    print('Making consensus of ' + str(len(comparelist2)) + ' sequences with similarity >= ' + str(similar))
    try:
        consensus = l_median(comparelist2)  # create consensuse sequence
        outfile = infile.split('.')[0] + '_consensus.fasta'
        with open(outfile, 'w')as of:
            of.write('>consensus\n' + consensus)
            print('Saving ' + outfile)
    except TypeError:
        print('!!!! No sequences have a similarity >= ' + str(similar) + ' !!!!')
#==============================================================================    
if __name__ == '__main__':
    args = get_arguments()
    infile = args.input  
    read_file(infile)