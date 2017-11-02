from difflib import SequenceMatcher
import unittest

class MgfItem:
  def __init__(self):
    self.Data = []

  def getPrecursorMz(self):
    for idx, val in enumerate(self.Data):
      if val.startswith("PEPMASS="):
        return float(val[8:])
    raise Exception("Cannot find PEPMASS!")
      
  def setPrecursorMz(self, precursorMz):
    for idx, val in enumerate(self.Data):
      if val.startswith("PEPMASS="):
        self.Data[idx] = "PEPMASS=%f" % precursorMz
        return
  
  def getTitle(self):
    for idx, val in enumerate(self.Data):
      if val.startswith("TITLE="):
        return val[6:]
    raise Exception("Cannot find TITLE!")
      
  def setTitle(self, title):
    for idx, val in enumerate(self.Data):
      if val.startswith("TITLE="):
        self.Data[idx] = "TITLE=%s" % title
        return
      
  def getCharge(self):
    for idx, val in enumerate(self.Data):
      if val.startswith("CHARGE="):
        return int(val[7])
    raise Exception("Cannot find CHARGE!")

  def read(self, sr):
    self.Data=[]
    for line in sr:
      if line.startswith("BEGIN IONS"):
        self.Data=[]
      elif line.startswith("END IONS"):
        return(True)
      else:
        self.Data.append(line.rstrip())
    return(False)
    
  def write(self, sw):
    sw.write("BEGIN IONS\n")
    for line in self.Data:
      sw.write("%s\n" % line)
    sw.write("END IONS\n\n")
    
class TestMgfItem(unittest.TestCase):
  def testTitle(self):
    item = MgfItem()
    item.Data = ["TITLE=B00_01_140613_zhuxu_Qu_HCD_OT_2hr.3.3.2.dta",
                 "PEPMASS=325.42654",
                 "CHARGE=2+",
                 "SCANS=3"]
    self.assertEqual("B00_01_140613_zhuxu_Qu_HCD_OT_2hr.3.3.2.dta", item.getTitle())
    
    item.setTitle("newTitle")
    self.assertEqual("newTitle", item.getTitle())

  def testCharge(self):
    item = MgfItem()
    item.Data = ["TITLE=B00_01_140613_zhuxu_Qu_HCD_OT_2hr.3.3.2.dta",
                 "PEPMASS=325.42654",
                 "CHARGE=2+",
                 "SCANS=3"]
    self.assertEqual(2, item.getCharge())

  def testPrecursorMz(self):
    item = MgfItem()
    item.Data = ["TITLE=B00_01_140613_zhuxu_Qu_HCD_OT_2hr.3.3.2.dta",
                 "PEPMASS=325.42654",
                 "CHARGE=2+",
                 "SCANS=3"]
    self.assertEqual(325.42654, item.getPrecursorMz())
    
    item.setPrecursorMz(300.05)
    self.assertEqual(300.05, item.getPrecursorMz())

