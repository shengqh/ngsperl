import sys

inputfile = sys.argv[1]
sequencefile = sys.argv[2]
outputfile = sys.argv[3]

#inputfile="/scratch/cqs/shengq1/vickers/20161121_smallRNA_3018_85_spikein_run2/human/class_independent/identical_sequence_count_table/pbs/human_spikein_sequence.filelist"
#sequencefile="/scratch/cqs/shengq1/vickers/20161121_smallRNA_3018_85_spikein_run2/spikein.txt"
#outputfile="/scratch/cqs/shengq1/vickers/20161121_smallRNA_3018_85_spikein_run2/human.txt"

dupcount_files=[]
with open(inputfile, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
          dupcount_files.append([parts[1], parts[0]])

sequences=[]
map={}
with open(sequencefile, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
          sequences.append([parts[0], parts[1]])
          map[parts[0]]={}

for dupfile in dupcount_files:
    name=dupfile[0]
    file=dupfile[1]
    with open(file, 'r') as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            seq=parts[2]
            for sequence in sequences:
                if sequence[1] == seq:
                    map[sequence[0]][name]=parts[1]

#print(map)

with open(outputfile, 'w') as f:
    f.write("Name\tSequence\t%s\n" % "\t".join([sample[0] for sample in dupcount_files]))
    for sequence in sequences:
        f.write("%s\t%s" % (sequence[0], sequence[1]))
        curmap = map[sequence[0]]
        for sample in dupcount_files:
            if sample[0] in  curmap:
                f.write("\t" + curmap[sample[0]])
            else:
                f.write("\t0")
        f.write("\n")

