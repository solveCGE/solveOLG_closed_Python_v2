import numpy as np

# formatted reporting
def report(reporttext, reportcalc):
  cursorstart = 45
  
  countlinebreaks = 0
  for i in range(len(reporttext)):
      if reporttext[:i+1] != "\n":
          break
      countlinebreaks += 1
  
  cursorstart = max(len(reporttext)-countlinebreaks,cursorstart)+3*(reportcalc>=0)+2*(reportcalc<0)
  
  charfill = " "*(cursorstart-len(reporttext)+countlinebreaks)
  
  print(reporttext+charfill+str(reportcalc))

def zeroscol(dim1):
  return np.zeros((dim1,1))

def zerosrow(dim1):
  return np.zeros(dim1)

def zerosmat(dim1,dim2):
  return np.zeros((dim1,dim2))

def onescol(dim1):
  return np.ones((dim1,1))

def onesrow(dim1):
  return np.ones(dim1)

def onesmat(dim1,dim2):
  return np.ones((dim1,dim2))

def coh2per(inmat):

  maxage, numcoh = inmat.shape
  numper = numcoh - (maxage-1)

  if numper <= 0:
    raise ValueError("coh2per: insufficient number of columns in input matrix")
  
  outmat = zerosmat(maxage, numper)
  
  for a in range(0, maxage):
    outmat[a, :] = inmat[a, (maxage-1-a):(numcoh-a)]
  
  return outmat

def aggcoh2per(inmat):
  return np.sum(coh2per(inmat), axis=0)

def per2coh(inmat):
  maxage, numper = inmat.shape
  numcoh = numper + (maxage-1)
  
  calibvec = inmat[:,0]
  
  outmat = zerosmat(maxage,numcoh)
  
  for a in range(0,maxage):
    if a < (maxage-1):
      outmat[a,0:maxage-1-a] = onesrow(maxage-a-1)*calibvec[a]
    outmat[a,maxage-1-a:numcoh-a] = inmat[a,:]
    if a > 0:
      outmat[a,numcoh-a:numcoh] = onesrow(a)*inmat[a,numper-1]
  
  return outmat


