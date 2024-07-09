print("Start calibration: \n")

# parameters
delta                    = 0.05                 # depreciation rate
r0                       = 0.04                 # real interest rate
sigma                    = 0.9                  # elasticity of inter-temporal substitution
sigL                     = 0.3                  # labor supply elasticity

# note: ages are off-set by 1 year, e.g. age group 1 contains 0-year olds
fag                      = 14                   # first economically active age-group (age 15)
rag0                     = 61.3                 # retirement age group (retirement age 62.3), non-whole numbers allowed
iag0                     = 51                   # first age group giving inter-vivo transfers

ivpc                     = 0.2                  # intervivo transfer received per capita

# some normalizations
N0                       = 100.0                # population
Y0                       = 100.0                # GDP
L0                       = 30.0                 # total labor supply in efficiency units
w0                       = 2.0                  # wage rate
Cons0                    = 0.55*Y0              # consumption share (calibrated using taul0)

### DEMOGRAPHY ###
for i in range(0,nag):
  gamv0[i] = 1-0.89**(nag-i)                    # some simple profile

# survival of last age group is 0
gamv0[nag-1] = 0
  
# compute demography
Nv0[0]  = 1
for i in range(1,nag):
  Nv0[i] = Nv0[i-1]*gamv0[i-1]

# rescale population
NB0             = 1/sum(Nv0)*N0
Nv0             = Nv0/sum(Nv0)*N0

avage0 = sum(Nv0.transpose()*arange(0,nag))/N0
report("REPORT: Average age:",avage0)
lifeexpect0 = lifeexpect(gamv0)
report("REPORT: Life-expectancy at age 0:", lifeexpect0[0,0])   ### HOW TO NOT NEEDING EXTRA INDEX??
report("REPORT: Life-expectancy at age 65:", lifeexpect0[65,0]) ### HOW TO NOT NEEDING EXTRA INDEX??

### AGE PROFILES ###

# indicator for not-retired
notretv0[0:floor(rag0)]   = 1                       # not retired
notretv0[floor(rag0)]     = rag0-floor(rag0)        # partly retired

# intervivo-transfers
ivv0[iag0-1:nag,0]            = -linspace(ivpc, ivpc*2, nag-iag0+1) # some increasing profile (from ivpc to 2*ivpc)
ivv0[fag-1:iag0-1]            = -sum(ivv0[iag0-1:nag]*Nv0[iag0-1:nag])/sum(Nv0[fag-1:iag0-1])*onescol(iag0-fag)

iv0                           = sum(ivv0*Nv0)
if abs(iv0)>1e-10:
  raise ValueError("ERROR: UNBALANCED INTERVIVO TRANSFERS!")
  
thetav0                     = zeroscol(nag)                               # labor productivity parameters
theta_peak                  = floor(rag0)-10                              # assumption: productivity peaks 10 years before retirement
thetav0[fag-1:theta_peak,0] = linspace(0.7,1,theta_peak-fag+1)
thetav0[theta_peak:nag,0]   = linspace(1,0.1,nag-theta_peak)

ellv0                       = L0/sum(Nv0*thetav0*notretv0)*onescol(nag)   # labor supply

# partition of population
Nc0   = sum(Nv0[0:fag-1])          # number of children
Nw0 	= sum(notretv0*Nv0)-Nc0 		 # number of workers
Nr0 	= sum((1-notretv0)*Nv0) 	   # number of retirees
report("REPORT: Old-age dependency ratio:",sum(Nv0[65:nag])/sum(Nv0[15:65]))
report("REPORT: Economic dependency ratio:",(Nc0+Nr0)/Nw0)
report("CHECK: Newborns - deaths:", sum((1-gamv0)*Nv0)-NB0)
report("CHECK: Children + workers + retriees - pop.:", Nc0+Nw0+Nr0-N0)

### POLICY PARAMETERS ###
tauWv0                   = 0.15*onescol(nag)            # wage tax rate worker & retiree
tauF0                    = 0.2                          # payroll tax rate
tauC0                    = 0.2                          # consumption tax rate
tauprof0                 = 0.1                          # profit tax rate
pv0                      = 0.65*sum(w0*ellv0*thetav0*Nv0)/N0*onescol(nag)  # old-age pension (65% of average wage earnings)
DG0                      = 0.6*Y0                       # government debt level (60% of GDP)

# cGv0 is used to balance budget in calibration
cGv0_profile             = 0.2*onescol(nag)
cGv0_profile[0:25,0]     = linspace(0.4,0.2,25)
cGv0_profile[54:nag,0]   = linspace(0.2,1.0,nag-55+1) # some U-shaped profile

# price of consumption and age specific prices and tax rates (but the same for all age groups)
pc0     = 1+tauC0
tauCv0  = tauC0*onescol(nag)
pcv0    = pc0*onescol(nag)
wv0     = w0*onescol(nag)
rv0     = r0*onescol(nag)

LS0     = sum(notretv0*ellv0*thetav0*Nv0) # aggregate labor supply
LD0     = LS0
uck0    = (r0+delta*(1-tauprof0))/(1-tauprof0) # user-cost of capital
K0      = (Y0-(1+tauF0)*w0*LD0)/uck0
Inv0    = delta*K0
alpha   = K0*uck0/(K0*uck0+LS0*((1+tauF0)*w0))
qTob0   = (1-tauprof0)*alpha*Y0/K0 + tauprof0*delta + (1-delta) # = 1+r0
TFP0    = Y0/((K0**alpha)*(LS0**(1-alpha)))
#LD0     = ((1-alpha)*TFP0/((1+tauF0)*w0))**(1/alpha)*K0 # also true
TaxF0   = tauprof0*(Y0-(1+tauF0)*w0*LD0-(delta*K0))
Div0    = Y0-(1+tauF0)*w0*LD0-Inv0-TaxF0
V0      = (1+r0)*Div0/r0

def calibfind(xcalib0):
  
  global rho, taul0, ab0, abv0, taulv0, cGv0, yv0, lambdav0, Consv0, Av0, Savv0
  global A0, P0, CG0, Exp0, tauW0, Rev0, PB0, edy0, edl0, eda0, edg0, ediv0, edab0
  
  retvar = np.zeros(5)
  
  rho      = xcalib0[0]
  cGscale  = xcalib0[1]
  taul0    = xcalib0[2]
  ab0      = xcalib0[3]
  lambdain = xcalib0[4]
  
  abv0[fag-1:nag]    = ab0/(N0-Nc0)*onescol(nag-fag+1) # children do not receive accidental bequest (workers start out with 0 assets)
  taulv0[fag-1:nag]  = taul0
  cGv0               = cGv0_profile+cGscale
  
  # INCOME
  yv0     = notretv0*(wv0*(1-tauWv0)*ellv0*thetav0)+(1-notretv0)*(1-tauWv0)*pv0-taulv0
  
  # CONSUMPTION FOR ALL AGE GROUPS
  
  # Euler equation
  lambdav0[fag-1] = lambdain
  for a in range(fag-1,nag-1):
    lambdav0[a+1] = lambdav0[a]/((1/(1+rho))*gamv0[a]*(1+rv0[a]))
  
  Consv0[fag-1:nag] = (pcv0[fag-1:nag]*lambdav0[fag-1:nag])**(-sigma)
  
  # assets
  Av0[fag-1] = 0
  for a in range(fag,nag):
     Av0[a]     = (1+rv0[a-1])*(Av0[a-1]+yv0[a-1]+ivv0[a-1]+abv0[a-1]-pcv0[a-1]*Consv0[a-1])
  
  Savv0   = Av0+yv0+ivv0+abv0-pcv0*Consv0
  
  # AGGREGATION
  A0      = sum(Av0*Nv0)                                # total assets
  P0      = sum((1-notretv0)*pv0*Nv0)                   # expend pensions
  CG0     = sum(cGv0*Nv0)                               # government consumption
  Exp0    = CG0+P0                                      # total primary expenditures
  tauW0   = sum(tauWv0*notretv0*ellv0*thetav0*Nv0)/LS0  # average wage tax rate
  Rev0    = TaxF0+(tauF0*LD0+tauW0*LS0)*w0+taul0*(Nw0+Nr0)+tauC0*Cons0+sum((1-notretv0)*tauWv0*pv0*Nv0) # total revenues
  PB0     = DG0*r0/(1+r0)                               # primary balance
  
  # EXCESS DEMANDS
  edy0    = Cons0+CG0+Inv0-Y0            # goods market
  edl0    = LD0-LS0                      # labor market
  eda0    = DG0+V0-A0                    # assets market
  edg0    = Rev0-Exp0-PB0                # government budget
  ediv0   = -iv0                         # intervivo transfers resource constraint
  edab0   = sum((1-gamv0)*Savv0*Nv0)-ab0 # accidental bequest resource constraint
    
  retvar[0]   = edy0
  retvar[1]   = edg0
  retvar[2]   = sum(Consv0*Nv0)-Cons0
  retvar[3]   = edab0
  retvar[4]   = Savv0[nag-1]
  
  return retvar

# MATCH CALIBRATION TARGETS;
xcalib0 = array([0.01, 0.3719, 0.40, 13, 1]) # starting guesses for fsolve()

xout = fsolve(calibfind,xcalib0,full_output=1,xtol=1e-10)
if xout[2] != 1:
  raise ValueError("NEWTON METHOD DID NOT CONVERGE!\n");

### CALIBRATION OF LABOR SUPPLY MARGINS ###

# set parl0 in order to reproduce ell0, FOC ell0
parlv0[fag-1:nag] = (wv0[fag-1:nag]*(1-tauWv0[fag-1:nag])*thetav0[fag-1:nag]/pcv0[fag-1:nag])*(ellv0[fag-1:nag]**(-1/sigL))*(Consv0[fag-1:nag]**(-1/sigma))
# set parl1 in order to normalize disutility of labor to 0
parlv1    = (sigL/(1+sigL))*parlv0*(ellv0**((1+sigL)/sigL))
dis_totv0 = (sigL/(1+sigL))*parlv0*(ellv0**((1+sigL)/sigL))-parlv1

report("REPORT: Asset-to-output ratio:", A0/Y0)

checkA0         = sum(Av0*Nv0)-A0
checkAv0        = Av0[nag-1]+yv0[nag-1]+ivv0[nag-1]+abv0[nag-1]-pc0*Consv0[nag-1] # end of period assets of last age group are zero
checkN0         = sum(Nv0)-N0

chkcalib = array([edy0,edl0,edg0,ediv0,eda0,edab0,checkA0,checkAv0[0],checkN0])

report("CHECK: Calibration:",sum(chkcalib))

# fill time-dependent variables with calibration values
Cons                          = Cons0*onesrow(tend)
DG                            = DG0*onesrow(tend)
Inv                           = Inv0*onesrow(tend)
LD                            = LD0*onesrow(tend)
LS                            = LS0*onesrow(tend)
K                             = K0*onesrow(tend)
N                             = N0*onesrow(tend)
NB                            = NB0*onesrow(tend)
PB                            = PB0*onesrow(tend)
TFP                           = TFP0*onesrow(tend)
ab                            = ab0*onesrow(tend)
pc                            = pc0*onesrow(tend)
r                             = r0*onesrow(tend)
rag                           = rag0*onesrow(tend)
tauC                          = tauC0*onesrow(tend)
tauF                          = tauF0*onesrow(tend)
tauW                          = tauW0*onesrow(tend)
taul                          = taul0*onesrow(tend)
tauprof                       = tauprof0*onesrow(tend)
uck                           = uck0*onesrow(tend)

# fill time-dependent and age-dependent variables with calibration values
Av                            = kron(Av0, onesrow(tend))
Az                            = kron(Av0, onesrow(ncoh))
Consv                         = kron(Consv0, onesrow(tend))
Consz                         = kron(Consv0, onesrow(ncoh))
Nv                            = kron(Nv0, onesrow(tend))
Nz                            = kron(Nv0, onesrow(ncoh))
Savv                          = kron(Savv0, onesrow(tend))
Savz                          = kron(Savv0, onesrow(ncoh))
abv                           = kron(abv0, onesrow(tend))
abz                           = kron(abv0, onesrow(ncoh))
cGv                           = kron(cGv0, onesrow(tend))
cGz                           = kron(cGv0, onesrow(ncoh))
ellv                          = kron(ellv0, onesrow(tend))
ellz                          = kron(ellv0, onesrow(ncoh))
gamv                          = kron(gamv0, onesrow(tend))
gamz                          = kron(gamv0, onesrow(ncoh))
ivv                           = kron(ivv0, onesrow(tend))
ivz                           = kron(ivv0, onesrow(ncoh))
lambdav                       = kron(lambdav0, onesrow(tend))
lambdaz                       = kron(lambdav0, onesrow(ncoh))
notretv                       = kron(notretv0, onesrow(tend))
notretz                       = kron(notretv0, onesrow(ncoh))
pv                            = kron(pv0, onesrow(tend))
pz                            = kron(pv0, onesrow(ncoh))
tauCv                         = kron(tauCv0, onesrow(tend))
tauCz                         = kron(tauCv0, onesrow(ncoh))
tauWv                         = kron(tauWv0, onesrow(tend))
tauWz                         = kron(tauWv0, onesrow(ncoh))
taulv                         = kron(taulv0, onesrow(tend))
taulz                         = kron(taulv0, onesrow(ncoh))
thetav                        = kron(thetav0, onesrow(tend))
thetaz                        = kron(thetav0, onesrow(ncoh))
rv                            = kron(rv0, onesrow(tend))
rz                            = kron(rv0, onesrow(ncoh))
wv                            = kron(wv0, onesrow(tend))
wz                            = kron(wv0, onesrow(ncoh))

if genplots:

  plt.figure()
  plt.plot(arange(0,nag),Av0)
  plt.xlabel("age")
  plt.ylabel("assets")
  
  plt.show()
  
  plt.figure()
  plt.plot(arange(0,nag),Consv0)
  plt.plot(arange(0,nag),notretv0*ellv0*thetav0*wv0*(1-tauWv0))
  plt.plot(arange(0,nag),(1-notretv0)*pv0*(1-tauWv0))
  plt.plot(arange(0,nag),cGv0)
  plt.xlabel("age")
  plt.legend(loc=2,labels=["consumption","net labor income","net pension income","public consumption"])

plt.show()
