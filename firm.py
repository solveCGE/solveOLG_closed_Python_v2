# production function and marginal products
def fY(K,L,TFP): return TFP*K**alpha*L**(1-alpha)
def MPK(K,L,TFP): return TFP*alpha*(K/L)**(alpha-1)
def MPKinv(MPKin,L,TFP): return L*(MPKin/(TFP*alpha))**(1/(alpha-1))
def MPL(K,L,TFP): return TFP*(1-alpha)*(K/L)**alpha

# invert the firm problem: finds r for given V (and LS, TFP and taxes)
def rdemand(assetsupply, maxiter = 20, tol = 1e-6, verbose = False):

  #import pdb pdb.set_trace()

  K2      = copy.copy(K)
  uck2    = copy.copy(uck)
  r2      = copy.copy(r)
  qTob2   = copy.copy(qTob)
  
  error   = inf
  iterr   = 0
  
  while True:
    
    iterr          += 1
    error_old      = error
    qTob_old       = qTob2
    
    K2             = assetsupply/qTob2
    #K2[0]          = K[0] #predetermined
    uck2           = MPK(K2,LD,TFP)
    qTob2          = (1-tauprof)*uck2 + tauprof*delta + (1-delta)
    
    error = sum(abs(qTob2-qTob_old))
    
    if verbose == True:
      print("Iteration:\t",iterr,"\t\tError:\t",error)
    
    if iterr > maxiter:
      if verbose:
        print("No convergence!!")
      break
    if error < tol:
      if verbose:
        print("Convergence!!")
      break
    if error > error_old:
      if verbose:
        print("Increasing error: stop at previous step!")
      qTob2 = qTob_old
      break

  r2[0:(tend-1)] = qTob2[1:tend]-1
  r2[tend-1] = r2[tend-2]
  
  return r2


