import os
import os.path
import logging
import errno
import shutil

def checkFileExists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

def readFileMap(fileName):
  checkFileExists(fileName)

  result = {}
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)

def checkFileMap(fileMap):
  for sname in fileMap.keys():
    sfile = fileMap[sname]
    checkFileExists(sfile)

def initializeLogger(logfile, name, loglevel=logging.DEBUG):
  logger = logging.getLogger(name)
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

def writeRScropt(outfilePrefix, rScript, optionMap={}):
  checkFileExists(rScript)

  targetScript = outfilePrefix + ".r"
  optionMap["outfilePrefix"] = outfilePrefix
  with open(targetScript, "wt") as fout:
    for key in optionMap.keys():
      fout.write("%s='%s'\n" % (key, optionMap[key]))

    fout.write("setwd('%s')\n" % os.path.dirname(os.path.abspath(targetScript)))
    
    with open(rScript, "rt") as fin:
      bFirstSetwd = True
      for line in fin:
        if line.startswith("setwd") and bFirstSetwd:
          bFirstSetwd = False
          continue

        bInOption = False
        for key in optionMap.keys():
          if line.startswith(key + "="):
            optionMap.pop(key)
            bInOption = True
            break

        if not bInOption:
          fout.write(line)

  return(targetScript)

def writeRmdScript(outfilePrefix, rmdScript, optionMap={}, copyRFunction=True, optionsToIndividualFile=True):
  checkFileExists(rmdScript)

  if copyRFunction:
    rFunScript = os.path.join( os.path.dirname(__file__), "Rfunctions.R")
    checkFileExists(rFunScript)

    targetFolder = os.path.dirname(os.path.abspath(outfilePrefix))
    targetRFunScript =  os.path.join(targetFolder, "Rfunctions.R")
    shutil.copyfile(rFunScript, targetRFunScript)

  if optionsToIndividualFile:
    optionToIndividualFileName = outfilePrefix + ".options"
    with open(optionToIndividualFileName, "wt") as fout:
      for key in sorted(optionMap.keys()):
        fout.write("%s\t%s\n" % (key, optionMap[key]))
    optionMap = {"option_file": os.path.basename(optionToIndividualFileName)}

  targetScript = outfilePrefix + ".rmd"
  with open(targetScript, "wt") as fout:
    with open(rmdScript, "rt") as fin:
      for line in fin:
        if line.startswith("```"):
          fout.write(line)
          for key in optionMap.keys():
            fout.write("%s='%s'\n" % (key, optionMap[key]))
          fout.write("\n")
          break
        else:
          fout.write(line)

      for line in fin:
        bInOption = False
        for key in optionMap.keys():
          if line.startswith(key + "=") or line.startswith(key + "<-"):
            optionMap.pop(key)
            bInOption = True
            break

        if not bInOption:
          fout.write(line)

  return(targetScript)
