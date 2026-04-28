import sys
import gzip
import os
import logging
import argparse
from tkinter import FALSE
import xml.etree.ElementTree as ET
import subprocess
import pandas as pd

DEBUG = False
NOT_DEBUG= not DEBUG

if DEBUG:
  genomeListFile="/scratch/stein_lab/shengq2/20200226_4233_4263_michelle_smallRNA_human_v5_byTiger/data_visualization/bacteria_count/result/StaRRA_human_4233_4263__fileList1.list"
  databaseFile = "/scratch/stein_lab/shengq2/20200226_4233_4263_michelle_smallRNA_human_v5_byTiger/nonhost_library/bowtie1_rRNA_pm_table/result/rRNA_pm_StaRRA_human_4233_4263.count.xml"
  taskReadFile = "/scratch/stein_lab/shengq2/20200226_4233_4263_michelle_smallRNA_human_v5_byTiger/data_visualization/reads_in_tasks/result/StaRRA_human_4233_4263.NonParallel.TaskReads.csv"
  outputFile="/scratch/stein_lab/shengq2/20200226_4233_4263_michelle_smallRNA_human_v5_byTiger/data_visualization/bacteria_count/result/StaRRA_human_4233_4263.tsv"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA count from count xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-g', '--genomeListFile', action='store', nargs='?', help='Input bacteria genome count xml list file', required=NOT_DEBUG)
  parser.add_argument('-d', '--databaseFile', action='store', nargs='?', help="Original rRNA database count xml file", required=FALSE)
  parser.add_argument('-t', '--taskReadFile', action='store', nargs='?', help="Task read count file", required=NOT_DEBUG)
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output count file", required=NOT_DEBUG)

  args = parser.parse_args()
  
  print(args)
  
  genomeListFile = args.genomeListFile
  databaseFile = args.databaseFile
  taskReadFile = args.taskReadFile
  outputFile = args.output

logger = logging.getLogger('getBacteriaCount')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def readFileList(fileName):
  df = pd.read_csv(fileName, sep='\t', header=None, usecols=[0], dtype=str)
  result = df[0].dropna().str.strip()
  result = result[result != ""]
  return result.tolist()
   
genomeFiles = readFileList(genomeListFile)

# Read all genome files into DataFrames and concatenate
genome_dfs = []
for genomeFile in genomeFiles:
  logger.info("Parsing " + genomeFile)
  df = pd.read_csv(genomeFile, sep='\t', index_col=0)
  genome_dfs.append(df)

if genome_dfs:
  # Verify all genome files have identical column names
  ref_columns = genome_dfs[0].columns.tolist()
  for i, df in enumerate(genome_dfs[1:], start=1):
    if df.columns.tolist() != ref_columns:
      raise ValueError(f"Column mismatch in genome file {genomeFiles[i]}: expected {ref_columns}, got {df.columns.tolist()}")
  result_df = pd.concat(genome_dfs)
else:
  result_df = pd.DataFrame()

# Parse XML for bacteria entries
if databaseFile:
  logger.info("Parsing " + databaseFile)
  xml_rows = []
  for event, elem in ET.iterparse(databaseFile, events=('end',)):
    if elem.tag != 'query':
      continue

    is_bacteria = False
    for loc in elem.findall('location'):
      if loc.get("seqname") == "Bacteria":
        is_bacteria = True
        break
    
    if is_bacteria:
      xml_rows.append((elem.get("seq"), elem.get("sample"), int(elem.get("count"))))

    elem.clear()

  if xml_rows:
    xml_df = pd.DataFrame(xml_rows, columns=["Sequence", "Sample", "Count"])
    logger.info("save non-host bacteria database data to db_pandas.txt for debugging")
    xml_df.to_csv("db_pandas.txt", sep='\t', index=False)
    xml_pivot = xml_df.pivot(index="Sequence", columns="Sample", values="Count").fillna(0).astype(int)
    # Combine genome and XML results
    result_df = result_df.combine_first(xml_pivot).fillna(0).astype(int)
  else:
    logger.info("No bacteria entries found in XML.")

# Remove duplicate sequences across files
result_df = result_df[~result_df.index.duplicated(keep='first')]

# Defragment after combine_first
result_df = result_df.copy()

# Sort by total count descending, then by sequence name alphabetically
result_df['_total'] = result_df.sum(axis=1)
result_df = result_df.reset_index()
result_df = result_df.sort_values(by=['_total', result_df.columns[0]], ascending=[False, True])
result_df = result_df.set_index(result_df.columns[0])
result_df = result_df.drop(columns='_total')

# Sort columns alphabetically
result_df = result_df.reindex(sorted(result_df.columns), axis=1)
result_df.index.name = "Sequence"

# Write main output
result_df.to_csv(outputFile, sep='\t')

# Write summary
summaryFile = outputFile + ".summary"
sample_totals = result_df.sum(axis=0)
summary_df = sample_totals.rename_axis("Sample").reset_index(name="Count")
summary_df.to_csv(summaryFile, sep='\t', index=False)

rscript = os.path.realpath(__file__) + ".R"
target_r = os.path.basename(rscript)
with open(target_r, "wt") as fout:
  fout.write("outFile='%s'\n" % summaryFile)
  fout.write("parFile1='%s'\n" % summaryFile)
  fout.write("parFile2='%s'\n" % taskReadFile)
  fout.write("setwd('%s')\n\n" % os.path.dirname(os.path.realpath(outputFile)))
  with open(rscript, "rt") as fin:
    for line in fin:
      fout.write(line)

subprocess.call("R --vanilla -f " + target_r, shell=True)

logger.info("done.")
