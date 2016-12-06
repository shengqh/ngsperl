import sys

inputfile = sys.argv[1]
outputfile = sys.argv[2]

#inputfile="Z:/Shared/Labs/Vickers Lab/ShilinZhao/20161003_smallRNA_3018-KCV-77_78_79_mouse_v2/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV-77_78_79.miRNA.NTA.count"
#outputfile="Z:/Shared/Labs/Vickers Lab/ShilinZhao/20161003_smallRNA_3018-KCV-77_78_79_mouse_v2/host_genome/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_KCV-77_78_79.miRNA.NTA.base.count"

samples=[]
ntacounts={}
with open(inputfile, 'r') as f:
    samples=[s.strip() for s in f.readline().split('\t')[3:]]
    for line in f:
        parts=line.strip().split('\t')
        if parts[0].endswith("NTA_"):
            continue
        ntas=parts[0].split("NTA_")[1]
        counts=[float(count_str) for count_str in parts[3:]]
        for i, c in enumerate(ntas):
            if c=='N':
               continue
            nta="%d%c"%(i+1,c)
            if nta in ntacounts:
                oldcount = ntacounts[nta]
                for i,c in enumerate(counts):
                    oldcount[i] = oldcount[i] + c
            else:
                ntacounts[nta]=list(counts)

with open(outputfile, 'w') as f:
    f.write("NTA\t%s\n" % "\t".join(samples))
    for nta in sorted(ntacounts):
        f.write("%s\t%s\n" % (nta, "\t".join([str(count) for count in ntacounts[nta]])))
