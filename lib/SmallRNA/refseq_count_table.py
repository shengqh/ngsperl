import logging
import re
import argparse
import copy
import gzip

def initialize_logger():
  logger = logging.getLogger('table')
  loglevel = logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

logger = initialize_logger()

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Get read count in refseq bacteria",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count list file', required=NOT_DEBUG)
parser.add_argument('-a', '--assembly', action='store', nargs='?', help='Input assembly summary file', required=NOT_DEBUG)
parser.add_argument('-t', '--taxonomy', action='store', nargs='?', help='Input species taxonomy file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output prefix", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/scratch/vickers_lab/projects/20220301_4893_2_RA_smRNA_mouse_v5_bacteria/count_table/result/bacteria__fileList1.list"
  args.assembly='/data/cqs/references/bacteria/20220302_assembly_summary_refseq.txt'
  args.taxonomy='/data/cqs/references/bacteria/20220224.taxonomy.txt'
  args.output="/scratch/vickers_lab/projects/20220301_4893_2_RA_smRNA_mouse_v5_bacteria/count_table/result/bacteria"

class Species:
  def __init__(self, name, taxonomy_id):
    self.taxonomy_id = taxonomy_id
    self.name = name
    self.queries = {}
    self.queries_set = {}
    self.is_subset = False
  
  def add_auery(self, sample, query, count):
    self.queries.setdefault(sample, {})[query] = count
    self.queries_set.setdefault(sample, set()).add(query)

  def sum_count(self):
    self.sample_count = {sample:sum(self.queries[sample].values()) for sample in self.queries.keys()}
    self.query_count = sum(self.sample_count.values())

  def merge(self, another):
    for sample in another.queries.keys():
      if sample not in self.queries:
        self.queries[sample] = {}
        self.queries_set[sample] = set()
      self.queries[sample].update(another.queries[sample])
      self.queries_set[sample].update(another.queries_set[sample])

  def contains(self, another):
    if self.query_count < another.query_count:
      return(False)
    
    for sample in another.sample_count.keys():
      if not sample in self.sample_count:
        return(False)
      if self.sample_count[sample] < another.sample_count[sample]:
        return(False)

    for sample in another.queries.keys():
      another_q = another.queries_set[sample]
      self_q = self.queries_set[sample]
      if not another_q.issubset(self_q):
          return(False)
      
    return(True)

logger.info("reading accession/taxonomy_id map from " + args.assembly + "...")
acc_taxoid_map = {}
with open(args.assembly, "rt") as fin:
  for line in fin:
    if line.startswith('#'):
      continue
    parts = line.split('\t')
    acc_taxoid_map[parts[0]] = parts[6]

logger.info("reading taxonomy/name map from " + args.assembly + "...")
taxoid_name_map = {}
with open(args.taxonomy, "rt") as fin:
  fin.readline()
  for line in fin:
    parts = line.split('\t')
    taxoid_name_map[parts[0]] = {
      'name':parts[2],
      'species':parts[11].rstrip(),
      'genus':parts[10],
      'family':parts[9],
      'order':parts[8],
      'class':parts[7],
      'phylum':parts[6],
      'kingdom':parts[5],
    }

count_map = {}
species_map = {}
samples=[]
with open(args.input, "rt") as fl:
  for line in fl:
    parts = line.rstrip().split('\t')
    count_file = parts[0]
    sample=parts[1]
    samples.append(sample)
    logger.info(f"parsing {count_file}")
    with gzip.open(count_file, "rt") as fin:
      fin.readline()
      for bl in fin:
        bparts = bl.split('\t')
        query = bparts[0].split(' ')[0]
        count = int(bparts[1])
        count_map[query] = count
        genbank_id_list = bparts[2].rstrip().split(',')
        for genbank_id in genbank_id_list:
          taxoid=acc_taxoid_map[genbank_id]
          name=taxoid_name_map[taxoid]['name']
          if name not in species_map:
            species_map[name] = Species(name, taxoid)
          species_map[name].add_auery(sample, query, count)

species_values = [v for v in species_map.values()]
for sv in species_values:
  sv.sum_count()

species_values.sort(key=lambda x:x.query_count, reverse=True)
# for i1 in range(0, len(species_values)-1):
#   if species_values[i1].is_subset:
#     continue
#   if i1 % 100 == 0:
#     logger.info(f"checking subset: {i1+1} / {len(species_values)}")
#   for i2 in range(i1+1, len(species_values)):
#     if species_values[i2].is_subset:
#       continue
#     if species_values[i1].contains(species_values[i2]):
#       species_values[i2].is_subset = True

# with open(args.output + ".species.txt", "wt") as fout:
#   fout.write("Feature\tTaxonomyId\t" + "\t".join(samples) + "\n")
#   for species in species_values:
#     countstr = "\t".join(str(species.sample_count[sample]) if sample in species.sample_count else "0" for sample in samples)
#     fout.write(f"{species.name}\t{species.taxonomy_id}\t{countstr}\n")

for level in [ 'species', 'genus', 'family', 'order', 'class', 'phylum' ]:
  logger.info(f"output {level} ...")
  level_map = {}
  for species in species_values:
    taxo_id = species.taxonomy_id
    level_id = taxoid_name_map[taxo_id][level]
    if level_id == "":
      level_id = '-1'
      level_name = 'Unclassfied'
    else:
      level_name = taxoid_name_map[level_id]['name']

    if level_id not in level_map:
      l = Species(level_name, level_id)
      l.queries = copy.deepcopy(species.queries)
      l.queries_set = copy.deepcopy(species.queries_set)
      level_map[level_id] = l
    else:
      l = level_map[level_id]
      l.merge(species)

  level_values = [v for v in level_map.values()]
  for sv in level_values:
    sv.sum_count()

  level_values.sort(key=lambda x:x.query_count, reverse=True)
  with open(f"{args.output}.{level}.txt", "wt") as fout:
    fout.write("Feature_name\tFeature_TaxonomyId\t" + "\t".join(samples) + "\n")
    for species in level_values:
      countstr = "\t".join(str(species.sample_count[sample]) if sample in species.sample_count else "0" for sample in samples)
      fout.write(f"{species.name}\t{species.taxonomy_id}\t{countstr}\n")

logger.info("done")
