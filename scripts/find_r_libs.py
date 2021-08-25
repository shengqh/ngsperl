import os
import os.path
from glob import glob

absfile=os.path.abspath(__file__)
PATH=os.path.dirname(os.path.dirname(absfile))
print(PATH)

rfiles=[y for x in os.walk(PATH) for y in glob(os.path.join(x[0], '*.r'))]
rfiles2=[y for x in os.walk(PATH) for y in glob(os.path.join(x[0], '*.R'))]
rmdfiles=[y for x in os.walk(PATH) for y in glob(os.path.join(x[0], '*.RMD'))]
rmdfiles2=[y for x in os.walk(PATH) for y in glob(os.path.join(x[0], '*.Rmd'))]
rmdfiles3=[y for x in os.walk(PATH) for y in glob(os.path.join(x[0], '*.rmd'))]

result=rfiles+rfiles2+rmdfiles+rmdfiles2+rmdfiles3

#print(result)
libs=set()
for rscript in result:
  with open(rscript, "rt") as fin:
    for line in fin:
      if "library(" in line:
        libs.add(line.strip())

alllibs=sorted([l for l in libs])
alllibs=[l.replace("library(", "") for l in alllibs]
alllibs=[l.replace(")", "") for l in alllibs]
alllibs=[l.replace('"', "") for l in alllibs]
alllibs=[l.replace("'", "") for l in alllibs]
alllibs=[l.replace(";", "") for l in alllibs]
alllibs=[l for l in alllibs if not l.startswith('#')]

alllibs.remove('mafreport')
alllibs.remove('ChIPQC')

alllibs.insert(0, "remotes")
alllibs.insert(0, "hdf5r")
alllibs.insert(0, "BiocManager")

github_libs = {
  "ChIPQC": 'shengqh',
  "cutoff": 'choisy',
  "mafreport": 'slzhao',
  "scRNABatchQC": 'liuqivandy'
}

final_libs = []
for l in alllibs:
  if l in github_libs:
    final_libs.append('%s/%s' % (github_libs[l], l))
  else:
    final_libs.append(l)

install_libs=[]
docker_libs=[]
for l in final_libs:
  if l == "BiocManager":
    install_str='install.packages'
  else:
    install_str='BiocManager::install'
  install_libs.append('if("' + l + '" %in% rownames(installed.packages()) == FALSE) {' + install_str + '("' + l + '", update=FALSE, ask=FALSE, dependencies=TRUE)}')
  docker_libs.append('RUN R -e "' + 'if(\'' + l + '\' %in% rownames(installed.packages()) == FALSE) {' + install_str + '(\'' + l + '\', update=FALSE, ask=FALSE, dependencies=TRUE)}"')

with open(os.path.dirname(absfile) + "/install_package.r", "wt") as fout:
  fout.write('Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")\n\n')
  for l in install_libs:
    fout.write(l + '\n')

with open(os.path.dirname(absfile) + "/install_package.docker", "wt") as fout:
  fout.write('RUN R -e "' + 'Sys.setenv(\'R_REMOTES_NO_ERRORS_FROM_WARNINGS\' = \'true\')"\n')
  for l in docker_libs:
    fout.write(l + '\n')
