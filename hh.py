# solve household problem for given prices wz[,z], abz[,z], taxes, etc.

@jit
def HH_root(lambdain, sage0, z0, lambdaz, pcz, Consz, ellz, dis_totz, yz, Az, Savz, gamz, rz, tauCz, wz, tauWz, thetaz, parlv0, parlv1, notretz, taulz, ivz, abz):

  sage = sage0 - 1
  z    = z0 - 1

  # EULER EQUATION: solve forward in age
  lambdaz[sage,z]    = lambdain
  if sage < nag:
    for a in range(sage,nag-1):
      lambdaz[a+1,z] = lambdaz[a,z]/((1/(1+rho))*gamz[a,z]*(1+rz[a,z]))

  # CONSUMPTION
  pcz[sage:nag,z]      = 1+tauCz[sage:nag,z]
  Consz[sage:nag,z]    = (pcz[sage:nag,z]*lambdaz[sage:nag,z])**(-sigma)

  # HOURS SUPPLY
  ellz[sage:nag,z]     = ((wz[sage:nag,z]*(1-tauWz[sage:nag,z])*thetaz[sage:nag,z]/pcz[sage:nag,z]*(Consz[sage:nag,z]**(-1/sigma)))/parlv0[sage:nag,0])**sigL
  dis_totz[sage:nag,z] = (sigL/(1+sigL))*parlv0[sage:nag,0]*ellz[sage:nag,z]**((1+sigL)/sigL)-parlv1[sage:nag,0]

  # CONSUMPTION AND SAVINGS
  yz[sage:nag,z]       = notretz[sage:nag,z]*(wz[sage:nag,z]*(1-tauWz[sage:nag,z])*ellz[sage:nag,z]*thetaz[sage:nag,z])+(1-notretz[sage:nag,z])*(1-tauWz[sage:nag,z])*pz[sage:nag,z]-taulz[sage:nag,z]

  # ASSETS: solve forward in age
  Az[0,z]         = 0
  
  if sage < nag:
    for a in range(sage,nag-1):
      Az[a+1,z]   = (1+rz[a,z])*(Az[a,z]+yz[a,z]+ivz[a,z]+abz[a,z]-pcz[a,z]*Consz[a,z]) # if sage > 1 take previous age entry in Az as starting value! (i.e. has to be given globally not passed in function)
  
  Savz[sage:nag,z]  = Az[sage:nag,z]+yz[sage:nag,z]+ivz[sage:nag,z]+abz[sage:nag,z]-pcz[sage:nag,z]*Consz[sage:nag,z]
  
  return Savz[nag-1,z]

def HH(sage0, z0, maxiter = 30, stol = 1e-10, atol = 0.1):

  err            = inf
  iterr          = 0
  trys           = 0
  stepsize       = 1e-6 # for numerical gradient
  
  lambdatrys     = array([1.0,0.5,1.5,0.25,1.25,0.1,1.0])
  maxtrys        = lambdatrys.size
  while_continue = True
  
  while while_continue:
    
    while_continue = False
    lambdazsave    = lambdaz[sage0-1,z0-1]
    
    while (err > stol or abs(Savz[nag-1,z0-1]) > atol) and trys < maxtrys:
      
      iterpertry = 0
      lambdaz1 = lambdazsave*lambdatrys[trys]; trys += 1
      
      breakwhile = False
      while err > stol and iterpertry < maxiter and breakwhile == False:
        if iterpertry == 0: # Newton step for first iteration
          f2 = HH_root(lambdaz1+stepsize, sage0, z0, lambdaz, pcz, Consz, ellz, dis_totz, yz, Az, Savz, gamz, rz, tauCz, wz, tauWz, thetaz, parlv0, parlv1, notretz, taulz, ivz, abz); iterr += 1
          
          if not isfinite(f2):
            breakwhile = True
            break
          f1 = HH_root(lambdaz1, sage0, z0, lambdaz, pcz, Consz, ellz, dis_totz, yz, Az, Savz, gamz, rz, tauCz, wz, tauWz, thetaz, parlv0, parlv1, notretz, taulz, ivz, abz); iterr +=  1
          
          if not isfinite(f1) or abs(f2-f1)<1e-16:
            breakwhile = True
            break
          lambdaz2 = lambdaz1 - f1*stepsize/(f2-f1)
          if not isfinite(lambdaz2) or lambdaz2 < 0:
            breakwhile = True
            break
        else: # Secant method
          f1 = HH_root(lambdaz1, sage0, z0, lambdaz, pcz, Consz, ellz, dis_totz, yz, Az, Savz, gamz, rz, tauCz, wz, tauWz, thetaz, parlv0, parlv1, notretz, taulz, ivz, abz); iterr += 1
          
          if not isfinite(f1) or abs(f1-f0)<1e-16:
            breakwhile = True
            break
          lambdaz2 = lambdaz1 - f1*(lambdaz1-lambdaz0)/(f1-f0)
          if not isfinite(lambdaz2) or lambdaz2 < 0:
            breakwhile = True
            break

        err = abs(lambdaz2-lambdaz1)
        lambdaz0   = lambdaz1
        lambdaz1   = lambdaz2
        f0         = f1
        iterpertry = iterpertry + 1

  if abs(Savz[nag-1,z0-1]) > atol:
    HH_nonconvz[nag-1,z0-1] = 1 # counter

def HHall(starttime = 1, calibinit = False, scaleA = 1):
  
  for z0 in range(starttime,ncoh+1):
    if z0 <= nag-fag+starttime-1:
      if calibinit == True:
        Az[:,z0-1] = Av0[:,0]
      Az[nag-1-(z0-starttime),z0-1] *= scaleA
      HH(nag-(z0-starttime), z0)
    else:
      HH(fag, z0)
