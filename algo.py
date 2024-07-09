# computes demographic transition (and updates intervivo transfers accordingly)
def compdemo():
  
  global Nv, Nz, N, Nc, ivv, ivz
  
  # compute demography transition
  for tt in range(1,tend):
    Nv[0,tt] = NB[tt]
    for i in range(1,nag):
      Nv[i,tt] = Nv[i-1,tt-1]*gamv[i-1,tt-1]

  Nz = per2coh(Nv)
  N  = aggcoh2per(Nz)
  Nc = sum(Nv[0:(fag-1),:],axis=0)
  
  # Compute neutral intervivo-transfers by rescaling received transfers
  for tt in range(0,tend):
    ivgiven          = -sum(Nv[:,tt]*ivv[:,tt]*(ivv[:,tt]<0))
    ivreceived       = sum(Nv[:,tt]*ivv[:,tt]*(ivv[:,tt]>0))
    
    ivv[ivv[:,tt]>0,tt] = ivv[ivv[:,tt]>0,tt]*(ivgiven/ivreceived)
    
    if abs(sum(ivv[:,tt]*Nv[:,tt]))>1e-10:
      print("ERROR IN RECOMPDEMO: Unbalanced intervivo transfers!")

  ivz = per2coh(ivv)


# computes the further life expectancy
def lifeexpect(gamv):
  
  nag = gamv.size
  
  lifeexpectageN = zeroscol(nag)
  for a in range(nag-2,-1,-1):
    lifeexpectageN[a] = (lifeexpectageN[a+1]+1)*gamv[a]+gamv[a+1]*(1-gamv[a])
  
  return lifeexpectageN

# main routine that solves the transition path of the full model
def solveOLG(starttime = 1, maxiter = 200, tol = 1e-4, damping_budget = 1.0, damping_assets = 1.0, damping_ab = 1.0, damping_r = 0.5, damping_new_assets = 0.7):
  
  global uck, K, Inv, qTob, Y, w, wz, V, TaxF, Div
  global Cons, LS, A, ab, iv, Nw, Nr, P, tauW, TaxP, Taxl, Rev, CG, Exp, PB
  global edy, edg, edl, eda, ediv, edab, edw
  global HH_nonconvz
  global tauWv, tauWz, tauF, tauC, tauCv, tauCz, taul, taulv, taulz, tauprof, cGv, CG
  global r, rz, rv, abv, abz, LD
  global Av, Consv, lambdav, Savv, dis_totv, ellv, ev, wv, pcv, yv
  
  print("\nRunning Tatonnement Algorithm for Transition:\n")
  
  tstart_loop       = default_timer()
  
  scaleA            = 1.0 # initialize
  scaleab           = onesrow(tend) # initialize

  #===== demography ======#
  compdemo() # recomputes demographic transition

  for iterr in range(1,maxiter+1):
    tstart_iter     = default_timer()

    #===== solve the firm problem for given labor demand ======#
    uck[starttime:tend]     = (r[starttime-1:tend-1]+delta*(1-tauprof[starttime:tend]))/(1-tauprof[starttime:tend])
    K[starttime:tend]       = MPKinv(uck[starttime:tend],LD[starttime:tend],TFP[starttime:tend])
    
    Inv[starttime-1:tend-1] = K[starttime:tend] - (1-delta)*K[starttime-1:tend-1]
    
    Inv[tend-1]     = delta*K[tend-1]
    qTob            = (1-tauprof)*MPK(K,LD,TFP) + tauprof*delta + (1-delta)
    
    Y               = fY(K,LD,TFP)
    w               = MPL(K,LD,TFP)/(1+tauF)
    wv              = kron(w,onescol(nag))
    wz              = per2coh(wv)
    V               = qTob*K
    TaxF            = tauprof*(Y-(1+tauF)*w*LD-delta*K)
    Div             = Y-(1+tauF)*w*LD-Inv-TaxF

    #===== solve the households' problem for given prices and tax rates ======#
    HHall(starttime = starttime, calibinit = (iterr == 1), scaleA = scaleA)

    #===== aggregation ======#
    Cons      = aggcoh2per(Consz*Nz)
    LS        = aggcoh2per(notretz*ellz*thetaz*Nz)
    A         = aggcoh2per(Az*Nz)
    ab        = aggcoh2per(abz*Nz)
    iv        = aggcoh2per(ivz*Nz) # should be 0 by construction
    Nw        = aggcoh2per(notretz*Nz)
    Nr        = aggcoh2per((1-notretz)*Nz)
    
    # government budget
    P           = aggcoh2per((1-notretz)*pz*Nz)
    tauW        = aggcoh2per(tauWz*notretz*ellz*thetaz*Nz)/LS
    TaxP        = aggcoh2per((1-notretz)*tauWz*pz*Nz)
    Taxl        = aggcoh2per(taulz*Nz)
    Rev         = TaxF+(tauF*LD+tauW*LS)*w+Taxl+tauC*Cons+TaxP
    CG          = sum(cGv*Nv,axis=0)
    Exp         = CG+P
    
    # follow given debt-path
    PB[starttime-1:tend-1]  = DG[starttime-1:tend-1]-DG[starttime:tend]/(1+r[starttime-1:tend-1])
    PB[tend-1]              = r[tend-1]*DG[tend-1]/(1+r[tend-1])
    
    #===== excess demands ======# 
    edy       = Inv+Cons+CG-Y
    edg       = Rev-Exp-PB
    edl       = LD-LS
    eda       = DG+V-A
    ediv      = -iv
    edab      = aggcoh2per((1-gamz)*Savz*Nz)-ab
    edw       = 1*edy + w*edl + ediv + edab + edg + eda - append(eda[1:tend],eda[tend-1])/(1+r) # Walras' Law

    # check Walras' Law: this always has to hold (even out of equilibrium)! If not there is something wrong with accounting in the model
    if max(abs(edw[starttime-1:tend-1]))> 1e-10:
      raise ValueError("Error: Walras Law does not hold!")
      
    tend_iter     = default_timer()

    #===== checking error and breaking loop ======# 	
    err             = sum(abs(edy[starttime-1:tend]))+sum(abs(edg[starttime-1:tend]))+sum(abs(edl[starttime-1:tend]))+sum(abs(eda[starttime-1:tend]))+sum(abs(ediv[starttime-1:tend]))+sum(abs(edab[starttime-1:tend]))
    err2            = log(err/tol)
    
    print("Iteration:  ", "{:>3}".format(iterr) ,"   scaleA: ", "{:.6f}".format(scaleA), "   scaleab: ", "{:.6f}".format(mean(scaleab)), "   non-conv.HH: ", "{:>2.0f}".format(sum(HH_nonconvz)), "   Time: ",  "{:.4f}".format(tend_iter-tstart_iter), " sec   log of err/tol: ", "{: .8f}".format(err2),sep="")
    
    if err2 < 0.0:
      print(" "*110 + "Convergence!\n")
      break
    if iterr == maxiter:
      print(" "*110 + "No Convergence!\n")
      break

    HH_nonconvz[:,:] = 0 # reset convergence counter
    
    #======= updating for next iteration =======#
    # budget rules
    budget_surplus  = edg*damping_budget;
    
    if budget_bal == 1:
      tauWv     = tauWv - kron(budget_surplus/(w*LS),onescol(nag))
      tauWz     = per2coh(tauWv)
    
    if budget_bal == 2:
      tauF       = tauF - budget_surplus/(w*LD) 
    
    if budget_bal == 3:
      tauC       = tauC - budget_surplus/Cons 
      tauCv      = kron(tauC,onescol(nag))
      tauCz      = per2coh(tauCv)
    
    if budget_bal == 4:
      taul               = taul - budget_surplus/(N-Nc)
      taulv[fag-1:nag,:] = kron(taul,onescol(nag-fag+1))
      taulz              = per2coh(taulv)
     
    if budget_bal == 5:
      tauprof    = tauprof - budget_surplus/(Y-(1+tauF)*w*LD-delta*K) 
    
    if budget_bal == 6:
      cGv        = cGv + kron(budget_surplus/N,onescol(nag))
      CG         = sum(cGv*Nv,axis=0)
    
    # price updating
    newassets       = damping_new_assets*(A-DG) + (1-damping_new_assets)*V
    r_new           = rdemand(newassets)
    r               = damping_r*r_new + (1-damping_r)*r
    rv              = kron(r,onescol(nag))
    rz              = per2coh(rv)

    scaleab         = 1+(aggcoh2per((1-gamz)*Savz*Nz)/ab-1)*damping_ab
    abv             = abv*kron(scaleab,onescol(nag))
    abz             = per2coh(abv)
    LD              = LS
    scaleA          = 1+((DG[starttime-1]+V[starttime-1])/A[starttime-1]-1)*damping_assets

  # convert cohort-view variables back to period-view variables
  # (those where only cohort-view variables were altered in solveOLG)
  Av       = coh2per(Az)
  Consv    = coh2per(Consz)
  lambdav  = coh2per(lambdaz)
  Savv     = coh2per(Savz)
  dis_totv = coh2per(dis_totz)
  ellv     = coh2per(ellz)
  pcv      = coh2per(pcz)
  yv       = coh2per(yz)
  
  tend_loop = default_timer()
  print("Computation time:\t", "{:.4f}".format(tend_loop-tstart_loop), "sec")
  print("CHECK SOLUTION:\t\t", max(abs(edy)+abs(edl)+abs(edg)+abs(eda)+abs(ediv)+abs(edab)))
  
  
