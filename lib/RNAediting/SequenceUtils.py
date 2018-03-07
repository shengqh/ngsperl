
def getMismatchIndex(seq1, seq2):
  minLength = min(len(seq1), len(seq2))
  result = [idx for idx in range(0, minLength) if seq1[idx] != seq2[idx]]
  return(result)

def calcIdenticalRate(seq1, seq2):
  minLength = min(len(seq1), len(seq2))
  numOfIdentical = minLength - len(getMismatchIndex(seq1, seq2))
  return(numOfIdentical * 1.0 / minLength)
