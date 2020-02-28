#!/usr/bin/env python
import collections
import csv
import os
import sys
import time
import operator

from optparse import OptionParser
from collections import defaultdict
from pyfaidx import Fasta

__author__ = "Sivasubramani S"
__email__ = "siva.subramani975@gmail.com"
__version__ = "0.0.1"

NEWLINE = '\n'
DELIM = '\t'
HEADERLINE = 'MarkerID' + DELIM \
    + 'Alleles' + DELIM \
    + 'countMajorAllele' + DELIM \
    + 'countMinorAllele' + DELIM \
    + 'countMissing'


"""
    This tool can help the user to filter the Hapmap file using many filtering criteira

"""

usage = "usage: python %prog \n\n\t" \
        "-i --inHMP\t- Input Hapmap file name" \
        "-o --outHMP\t- Output Hapmap file name"

parser = OptionParser(usage=usage)
parser.add_option("-i", "--inHMP", dest="inHMP",
                  help="Input hapmap file",
                  metavar="<FILENAME>")
parser.add_option("-o", "--outFile", dest="outFile", help="Output Summary file",
                  metavar="<string>")
parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()


def initialize_2nucl():
    """
    initializing iupac nucleotide dictionary
    :return:
    """
    bases = defaultdict(dict)
    bases['A']['A'] = 'A'
    bases['T']['T'] = 'T'
    bases['G']['G'] = 'G'
    bases['C']['C'] = 'C'
    bases['N']['N'] = 'N'
    bases['A']['T'] = 'W'
    bases['T']['A'] = 'W'
    bases['A']['G'] = 'R'
    bases['G']['A'] = 'R'
    bases['A']['C'] = 'M'
    bases['C']['A'] = 'M'
    bases['T']['G'] = 'K'
    bases['G']['T'] = 'K'
    bases['T']['C'] = 'Y'
    bases['C']['T'] = 'Y'
    bases['G']['C'] = 'S'
    bases['C']['G'] = 'S'
    return bases


def initialize_1nucl():
    """
    initializing iupac 2 letter nucleotide dictionary
    :return:
    """
    bases = defaultdict(dict)
    bases['A'] = ['A', 'A']
    bases['T'] = ['T', 'T']
    bases['G'] = ['G', 'G']
    bases['C'] = ['C', 'C']
    bases['N'] = ['N', 'N']
    bases['W'] = ['A', 'T']
    bases['R'] = ['A', 'G']
    bases['M'] = ['A', 'C']
    bases['K'] = ['G', 'T']
    bases['Y'] = ['C', 'T']
    bases['S'] = ['C', 'G']
    return bases


def readColumns(row):
    """
    returning data columns from hapmap entry
    :param row:
    :return:
    """
    genotypes = row[11:]
    return genotypes


"""__main__"""

required = "inHMP outFile".split()
for req in required:
    if options.__dict__[req] is None:
        parser.error("Required option %s missing" % req)
inHMP = options.inHMP
outFile = options.outFile
lineNo = 0

startTime = int(time.time())

iupac2nucl = initialize_1nucl()


def initialize_baseCount():
    """
    initialize baseCount dictionary for each hapmap entry
    :return:
    """
    baseCount = defaultdict(dict)
    baseCount['A'] = 0
    baseCount['T'] = 0
    baseCount['G'] = 0
    baseCount['C'] = 0
    baseCount['N'] = 0
    return baseCount


def sortDictByValue(inDict, ascend):
    """
    Sort the dictionary based on its values
    :param inDict:
    :param ascend:
    :return:
    """
    sortedDict = defaultdict(dict)
    if ascend:
        for w in sorted(inDict, key=inDict.get, reverse=True):
            sortedDict[w] = inDict[w]
    else:
        for w in sorted(inDict, key=inDict.get, reverse=False):
            sortedDict[w] = inDict[w]
    return sortedDict


def getBaseCount(genotypeCalls):
    """
    get number of alleles from Hapmap data entry
    :param genotypeCalls:
    :return:
    """
    countBase = initialize_baseCount()
    for allele in genotypeCalls:
        for base in iupac2nucl[allele]:
            countBase[base] += 1
    return countBase


with open(inHMP) as iHandle:
    outFileHandle = open(outFile, 'w')
    for row in iHandle:
        row = row.rstrip()
        rowEntries = row.split(DELIM)
        lineNo += 1
        if lineNo == 1:
            genotypeNames = readColumns(rowEntries)
            outFileHandle.write(HEADERLINE)
            outFileHandle.write(NEWLINE)
        else:
            genotypeCalls = readColumns(rowEntries)
            markerId = rowEntries[0]
            majorAllele = ''
            minorAllele = ''
            thirdAllele = ''
            countBase = getBaseCount(genotypeCalls)
            countBase = sortDictByValue(countBase, True)
            i = 0
            for base, count in countBase.items():
                if base == 'N':
                    continue
                if count == 0:
                    continue
                else:
                    i += 1
                    if i == 1:
                        majorAllele = base
                    if i == 2:
                        minorAllele = base
                    if i == 3:
                        thirdAllele = base
            if majorAllele == '':
                print(markerId + " contains only missing calls. Skipping")
                continue
            if minorAllele == '':
                print(markerId + " is a monomorphic marker. Skipping")
                continue
            if thirdAllele != '':
                print(markerId + " has more then 2 alleles. Skipping")
                continue
            if countBase.get(majorAllele) == countBase.get(minorAllele):
                tmpMajor = min(majorAllele, minorAllele)
                tmpMinor = max(majorAllele, minorAllele)
                minorAllele = tmpMinor
                majorAllele = tmpMajor
            outFileHandle.write(markerId + DELIM \
                                + majorAllele + "/" + minorAllele + DELIM \
                                + str(countBase.get(majorAllele)) + DELIM \
                                + str(countBase.get(minorAllele)) + DELIM \
                                + str(countBase.get('N')))
            outFileHandle.write(NEWLINE)
    print(lineNo)
outFileHandle.close()
endTime = int(time.time())

print(endTime - startTime)
