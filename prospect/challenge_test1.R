# Libraries used
library(ProSpect)
library(LaplacesDemon)
library(foreach)
library(celestial)
library(magicaxis)

# Data used
data_path       = './Projects/Challenge/Data/'
stripe82_file   = 'Stripe82_highs2n_10subsample.csv'
stripe82_ssdata = read.csv(file.path(data_path, stripe82_file))
stripe82_df     = data.frame(stripe82_ssdata)[1,]              # testing 1 obj

print(stripe82_df$photoZ)

print(stripe82_df)

# Load the stellar populations template
bc03 <- data("BC03lr")
vz16 <- data("EMILES")

# The filters we will be using
filters = c('u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'z_SDSS')

filtout=foreach(i = filters)%do%{approxfun(getfilt(i))}

temppiv=pivwave[pivwave$filter %in% filters,]


# Input star-formation history (SFH) parameters
inpar=c(mSFR = 0,            #log-space
        mpeak = 0.7,         #log-space
        mperiod = 0.3,       #log-space
        mskew = 0.3,
        tau_birth = 0,       #log-space
        tau_screen = -0.5,   #log-space
        alpha_SF_birth = 1,
        alpha_SF_screen = 3
)

plot_sfh = function(par,agemax=13.3, add=FALSE, col='black', ylim=NULL,...){
  magcurve(massfunc_snorm_trunc(age=x, mSFR=10^par[1], mpeak=10^par[2], 
                                mperiod=10^par[3], mskew=par[4], 
                                magemax=agemax), 0, 13.8e9, add=add, col=col, 
           ylim=ylim,xlab='Age (Yr)', ylab='SFR (Msol / Yr)',...)
}

plot_sfh(inpar)

# 

gen_sed=ProSpectSED(massfunc=massfunc_snorm_trunc,
                   mSFR=10^inpar[1],
                   mpeak=10^inpar[2],
                   mperiod=10^inpar[3],
                   mskew=inpar[4],
                   tau_birth=10^inpar[5], 
                   tau_screen=10^inpar[6], 
                   alpha_SF_birth=inpar[7], 
                   alpha_SF_screen=inpar[8],
                   z=stripe82_df$photoZ,
                   Z=Zfunc_massmap_lin,
                   filtout=filtout,
                   Dale=Dale_NormTot,
                   speclib=BC03lr,
                   agemax=agemax
)
