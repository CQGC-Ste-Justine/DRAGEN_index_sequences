#!/bin/env python

IndexBarcodes = '/staging2/data/index_sequences/sequencing_index.processed.ALL.txt'

import sys, getopt, csv, os

def main(argv):
  inputdir = ''
  outputfile = ''
  try:
     opts, args = getopt.getopt(argv,"hd:o:",["idir=","ofile="])
  except getopt.GetoptError:
    print 'index_stats.py -d <Reports_dir>'
    sys.exit(42)
  for opt, arg in opts:
    if opt == '-h':
      print 'index_stats.py -d <Reports_dir>'
      sys.exit()
    elif opt in ("-d", "--idir"):
      inputdir = arg
    elif opt in ("-o", "--ofile"):
      outputfile = arg
      print "outputfile not implemented yet!"
      sys.exit(42)
#  print 'Input dir is', inputdir
#  print 'Output file is', outputfile
  if inputdir == '': 
    print ("Invalid inputdir.")
    print 'index_stats.py -d <Reports_dir>'
    sys.exit(42)
# ================================================== #

  scriptdir=os.path.abspath(os.path.dirname(sys.argv[0]))
  IndexBarcodes = scriptdir + '/sequencing_index.processed.ALL.txt'
  index1barcodes = {}
  index2barcodes = {}
  indexBarcodes = {}
  #BarcodeName,index1,index2
  #KAPA_V1_UDI1_A01,GTAACATC,AATCGCTG
  with open(IndexBarcodes, "r") as IndexBarcodesFH:
    for line in IndexBarcodesFH:
      BarcodeName, index1, index2 = line.strip().split(',')
      index1barcodes[index1] = BarcodeName
      index2barcodes[index2] = BarcodeName
      indexBarcodes[(index1,index2)] = BarcodeName
#  print ("indexBarcodes: " + str(indexBarcodes))

# ================================================== #

  SampleBarcode={}
  BarcodeCount={}
  UndeterminedCount=0
  #Lane,SampleID,Index,# Reads,# Perfect Index Reads,# One Mismatch Index Reads,# of >= Q30 Bases (PF),Mean Quality Score (PF)
  #1,15573,CTACCAGG-ATCAGTTG,6104339,5956580,147759,877964362,34.53
  with open (inputdir+"/Demultiplex_Stats.csv", "r") as DemultiplexStatsFileFH:
    next(DemultiplexStatsFileFH)
    for line in DemultiplexStatsFileFH:
      Lane, SampleID, Index, Reads, PerfectIndexReads,NumOneMismatchIndexReads,NumQ30Bases,MeanQualityScorePF = line.strip().split(',')
      if SampleID != 'Undetermined':
        index1, index2 = Index.split('-')
        indexes = (index1,index2)
        if indexBarcodes.has_key(indexes):
          barcodeName = indexBarcodes[indexes]
          SampleBarcode[barcodeName] = SampleID
          if BarcodeCount.has_key(barcodeName):
              BarcodeCount[barcodeName] += int(Reads)
          else: BarcodeCount[barcodeName] = int(Reads)
      else: 
        UndeterminedCount += int(Reads)
#  print ("SampleBarcode: " + str(SampleBarcode))
#  print ("BarcodeCount: " + str(BarcodeCount))

# ================================================== #

  #Lane,index,index2,# Reads
  #1,ACACGATC,ATGGTATT,6434455
  with open (inputdir+"/Top_Unknown_Barcodes.csv", "r") as TopUnknownBarcodesFH:
    next(TopUnknownBarcodesFH)
    for line in TopUnknownBarcodesFH:
      Lane, index1, index2, Reads = line.strip().split(',')
      indexes = (index1,index2)
      if indexBarcodes.has_key(indexes):
        barcodeName = indexBarcodes[indexes]
        UndeterminedCount -= int(Reads)
        if BarcodeCount.has_key(barcodeName):
          BarcodeCount[barcodeName] += int(Reads)
        else: BarcodeCount[barcodeName] = int(Reads)

  BarcodeCount['Unknown'] = UndeterminedCount
  SampleBarcode['Unknown'] = 'Unknown'
#  print ("SampleBarcode: " + str(SampleBarcode))
#  print ("BarcodeCount: " + str(BarcodeCount))

# ================================================== #

  writer = csv.writer(sys.stdout)
  writer.writerow(['#Sample', 'Barcode', 'Count', 'PercentTotal'])
  
  totalReads=sum(BarcodeCount.values())

#  for barcodeName in BarcodeCount.keys():
  for barcodeName in sorted(BarcodeCount, key=BarcodeCount.get, reverse=True):
    if SampleBarcode.has_key(barcodeName):
      sampleName = SampleBarcode[barcodeName]
    else: sampleName = 'Unknown'
    percent = float(BarcodeCount[barcodeName])/float(totalReads) * 100
    writer.writerow([sampleName, barcodeName, BarcodeCount[barcodeName], '%.2f' % percent])

# ================================================== #

if __name__ == "__main__":
  main(sys.argv[1:])
