# testing the examples here: 
# https://github.com/asgr/ProSpect/blob/master/vignettes/QuickFit.Rmd

library(ProSpect)
library(LaplacesDemon)
library(foreach)
library(celestial)
library(magicaxis)

bc03 <- data("BC03lr")
dale <- data("Dale_NormTot")
wl_pivot <- data("pivwave")

set.seed(0)
redshift = 0.1
filters = c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 
            'z_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 
            'W1_WISE' , 'W2_WISE', 'W3_WISE', 'W4_WISE', 'P100_Herschel', 
            'P160_Herschel', 'S250_Herschel' , 'S350_Herschel', 'S500_Herschel')

filtout=foreach(i = filters)%do%{approxfun(getfilt(i))}

temppiv=pivwave[pivwave$filter %in% filters,]

agemax=13.3e9-cosdistTravelTime(z=redshift, H0 = 67.8, OmegaM = 0.308)*1e9

inpar=c(mSFR = 0,            #log-space
        mpeak = 0.7,         #log-space
        mperiod = 0.3,       #log-space
        mskew = 0.3,
        tau_birth = 0,       #log-space
        tau_screen = -0.5,   #log-space
        alpha_SF_birth = 1,
        alpha_SF_screen = 3
)

plotSFH=function(par,agemax=13.3, add=FALSE,col='black',ylim=NULL,...){
  magcurve(massfunc_snorm_trunc(age=x,mSFR=10^par[1],mpeak=10^par[2],mperiod=10^par[3],
                                mskew=par[4], magemax=agemax),0,13.8e9,add=add,col=col,
           ylim=ylim,xlab='Age (Yr)', ylab='SFR (Msol / Yr)',...)
}

plotSFH(inpar)


genSED=ProSpectSED(massfunc=massfunc_snorm_trunc,
                   mSFR=10^inpar[1],
                   mpeak=10^inpar[2],
                   mperiod=10^inpar[3],
                   mskew=inpar[4],
                   tau_birth=10^inpar[5], 
                   tau_screen=10^inpar[6], 
                   alpha_SF_birth=inpar[7], 
                   alpha_SF_screen=inpar[8],
                   z=0.1,
                   Z=Zfunc_massmap_lin,
                   filtout=filtout,
                   Dale=Dale_NormTot,
                   speclib=BC03lr,
                   agemax=agemax
)


flux_input=data.frame(filter=temppiv$filter, pivwave=temppiv$pivwave, flux=genSED$Photom, fluxerr=genSED$Photom*0.1)
print(flux_input)

LumDist_Mpc = cosdistLumDist(z=0.1, H0 = 67.8, OmegaM = 0.308)


Data=list(flux=flux_input,
          arglist=list(z=0.1, massfunc=massfunc_snorm_trunc, agemax=agemax, Z=Zfunc_massmap_lin, LumDist_Mpc=LumDist_Mpc),
          speclib=BC03lr, 
          Dale=Dale_NormTot, 
          filtout=filtout, 
          SFH=SFHfunc, # the preferred functional form of the SFH (eg either SFHfunc, SFHburst)
          parm.names=c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen',
                       'alpha_SF_birth','alpha_SF_screen'), # which parameters to fit for
          logged=c(T,T,T,F,T,T,F,F), # fit parameters in logged or linear space
          intervals=list(lo=c(-4,-2,-1,-0.5,-2.5,-2.5,0,0), hi=c(3,1,1,1,1.5,1,4,4)), # fitting range for parameters
          fit = 'LD', # specifies the way in which the SED should be fitted ('LD', 'optim', 'CMA', or 'check')
          mon.names=c('LP','masstot','SFRburst',paste('flux.',flux_input$filter,sep='')),
          N=length(filters), # number of observed filters
          like='norm',
          verbose=FALSE
)

set.seed(1)
LDout = LaplacesDemon(Model=ProSpectSEDlike, Data=Data,  Initial.Values=inpar,
                      control=list(abstol=0.1), Iterations=1e4, Algorithm='CHARM', Thinning=1)

set.seed(1)
Data$fit = 'optim'
optout = optim(par=inpar, fn=ProSpectSEDlike, Data=Data, method = 'BFGS')



set.seed(1)
library(cmaeshpc)
Data$fit = 'CMA'
badpar = (Data$intervals$lo + Data$intervals$hi) / 2 #CMA is pretty tolerant of terrible initial guesses, unlike optim and LD.
CMAout = cmaeshpc(par=badpar, fn=ProSpectSEDlike, Data=Data, lower=Data$intervals$lo,
                  upper=Data$intervals$hi, control=list(trace=TRUE, maxwalltime=2))
print(CMAout$par)

maghist(LDout$Monitor[,"masstot"], verbose = FALSE, xlab='Stellar Mass / Msol', ylab='PDF')
abline(v=genSED$Stars$masstot, col='red')


magplot(flux_input$pivwave, LDout$Monitor[1,4:23], type='l', log='xy', grid=TRUE,
        xlab="Wavelength (Ang)", ylab='Flux Density / Jy')
for(i in 2:1e4){
  lines(flux_input$pivwave, LDout$Monitor[i,4:23], col=hsv(alpha=0.1))
}                                                       # all the 10^4 results!! 
points(flux_input[,c("pivwave","flux")])
magerr(flux_input$pivwave, flux_input$flux, ylo=flux_input$fluxerr)

colvec=hsv(h=magmap(BC03lr$Age,lo=1e6,hi=1e10,range=c(0,2/3), stretch = 'log', flip=T, type='num', bad=2/3)$map, v=0.7, alpha=0.5)


for(j in 1:6){
  par(mar=c(3.1,3.1,1.1,1.1))
  magplot(BC03lr$Wave, BC03lr$Zspec[[j]][1,], type='l', log='xy', col=colvec[1], xlim=c(3e2,1e5), ylim=c(1e-7,10), majorn=c(5,8), xlab=BC03lr$Labels$Wavelab, ylab=BC03lr$Labels$Lumlab, grid=TRUE)
  for(i in 2:221){lines(BC03lr$Wave, BC03lr$Zspec[[j]][i,], col=colvec[i])}
  legend('topleft', legend=paste('Z=',BC03lr$Z[j]))
  magbar('topright',range = c(1e6,1e10),log = TRUE, col = hsv(h = seq(2/3, 0, len = 10), v=0.7), title = 'Age / Yrs', titleshift = 0.5)
}
