import sys
import gzip
import os
import logging
import argparse
import string
import subprocess
import multiprocessing
from coltron import utils

# Set up parallel computing
num_cores = multiprocessing.cpu_count()
subprocess.call('ulimit -s unlimited', shell=True)

# bamToGFF_turbo.py

'''
The MIT License (MIT)

Copyright (c) 2013 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

# 20130716

# script to grab reads from a bam that align to a .gff file

# uses the bamliquidator super fast uber thingy written by Xin Zhou

def mapBamToGFF(bamliquidatorLocation, bamFile, gff, sense='.', extension=200, rpm=False, nBin=200):
  '''maps reads from a bam to a gff'''

  # creating a new gff to output
  newGFF = []
  # reading in the bam
  bam = utils.Bam(bamFile)

  # getting RPM normalization
  if rpm:
      MMR = round(float(bam.getTotalReads('mapped')) / 1000000, 4)
  else:
      MMR = 1

  print('using a MMR value of %s' % (MMR))

  # creating a sense trans
  senseTrans = string.maketrans('-+.', '+-+')

  # reading in the gff
  if type(gff) == str:
      gff = utils.parseTable(gff, '\t')

  # setting up a maxtrix table
  bamName = bamFile.split('/')[-1]
  newGFF.append(['GENE_ID', 'locusLine'] + ['bin_' + str(n) + '_' + bamName for n in range(1, nBin + 1, 1)])

  # getting and processing reads for gff lines
  ticker = 0
  print('Number lines processed')
  for line in gff:
    line = line[0:9]
    if ticker % 100 == 0:
        print(ticker)
    ticker += 1
    gffLocus = utils.Locus(line[0], int(line[3]), int(line[4]), line[6], line[1])

    # get the nBin and binSize
    binSize = gffLocus.len() / nBin
    # some regions will be too short to get info on
    if binSize == 0:
        clusterLine = [gffLocus.ID(), gffLocus.__str__()] + ['NA'] * nBin
        newGFF.append(clusterLine)
        continue

    # flippy flip if sense is negative
    if sense == '-':
        bamSense = string.translate(gffLocus.sense(), senseTrans)
    elif sense == '+':
        bamSense = gffLocus.sense()
    else:
        bamSense = '.'
    # using the bamLiquidator to get the readstring
    # print('using nBin of %s' % nBin)
    bamCommand = "%s %s %s %s %s %s %s %s" % (bamliquidatorLocation, bamFile, line[0], gffLocus.start(), gffLocus.end(), bamSense, nBin, extension)
    # print(bamCommand)
    getReads = subprocess.Popen(bamCommand, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    readString, stderr = getReads.communicate()
    if stderr:
        print("STDERR out: %s" % (stderr))
    denList = readString.split('\n')[:-1]
    # print("denlist is: %s" % denList)
    # flip the denList if the actual gff region is -
    if gffLocus.sense() == '-':
        denList = denList[::-1]

    # converting from units of total bp of read sequence per bin to rpm/bp

    denList = [round(float(x) / binSize / MMR, 4) for x in denList]

    # if the gff region is - strand, flip the

    clusterLine = [gffLocus.ID(), gffLocus.__str__()] + denList
    newGFF.append(clusterLine)

  return newGFF

def getName(fileName):
    return fileName.split('/')[-1].split('.')[0]
  
def mapBams(bamliquidatorLocation, dataDict, outputPrefix, gff, extension=200, nBin=200, overWrite=False, rpm=True):
  nameList = dataDict.keys()

  listFile = outputPrefix + ".gff.filelist"
  with open(listFile, "w") as fw:
    for name in nameList:
      print ('mapping %s to gff' % name)
      bamFile = dataDict[name]
      outFile = outputPrefix + '.' + name + '.gff'
      if os.path.isfile(outFile) and not overWrite:
        print ("File exists, ignored :" + outFile)
        continue
      
      newGFF = mapBamToGFF(bamliquidatorLocation, bamFile, gff, '.', extension, rpm, nBin)
      utils.unParseTable(newGFF, outFile, '\t')
      fw.write("%s\t%s\n" % (outFile, name))

def main():
  appPath = os.path.dirname(os.path.abspath(__file__))
  sys.path.append(appPath)
  
  from bedCenter import centerPeaks
  
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Generate heatmap for bam files based on assigned peak list.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('--bamliquidator', action='store', nargs='?', required=NOT_DEBUG, help='bamliquidator location')
  parser.add_argument('-i', '--input', action='store', nargs='?', required=NOT_DEBUG, help="Input peak file in bed format")
  parser.add_argument('-d', '--dataFile', action='store', nargs='?', required=NOT_DEBUG, help="Data file contains bam information without header, first column is bam file and second column is bam name")
  parser.add_argument('-o', '--output', action='store', nargs='?', required=NOT_DEBUG, help="Output prefix")
  parser.add_argument('-w', '--window', action='store', nargs='?', required=NOT_DEBUG, help="Peak window in each side")
  parser.add_argument('-e', '--extension', action='store', nargs='?', required=NOT_DEBUG, help="Read extension")

  args = parser.parse_args()
  
  if(DEBUG):
    args.bamliquidator = "/scratch/cqs/shengq1/local/bin/bamliquidator"
    args.input = "/scratch/cqs/shengq1/brown/20170327_chipseq_3593_25to30_human/macs1callpeak/result/No_Tx_H3K27ac/No_Tx_H3K27ac_peaks.bed"
    args.dataFile = "/scratch/cqs/shengq1/brown/heatmap/bam.filelist"
    args.output = "/scratch/cqs/shengq1/brown/heatmap/test"
    args.window = 5000
    args.extension = 200
  
  window=int(args.window)
  
  print(args)
  
  dataDict = {}
  with open(args.dataFile, "r") as f:
    for line in f:
      parts = line.rstrip().split('\t')
      if(len(parts) > 1):
        dataDict[parts[1]] = parts[0]

  centerBed = args.output + ".center.bed"
  centerPeaks(args.input, window, centerBed)

  bed = utils.parseTable(centerBed, '\t')
  gff = utils.bedToGFF(bed)

  mapBams(args.bamliquidator, dataDict, args.output, gff, extension=args.extension, nBin=200, overWrite=False, rpm=True)

if __name__ == "__main__":
    main()
