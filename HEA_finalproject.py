## Here is some python code just to demonstrate a few techniques you may need to use in order to write your code. This is not a working code in itself (it has no overall structure), but it just gives some examples of kinds of commands/routines you may need to use when you write your own code.
    
# Here are some packages you may be likely to use
    
## Import relevant packages:

from numpy import array, float
from random import random
from math import pow, sqrt, pi, log, sin, cos
from pylab import *


# First, in Python all you need to call a random fraction from 0 to 1 is the command random(). You can then multiply this by any limits to generate a random distribution of any variable. Eg.:
    
    theta_rand = random() * 2*pi
    
    # That is simply a random angle between 0 and 2*pi. So any time you need to select a value from some kind of distribution, all you need to do is multiply a randomly generated fraction (via the 'random()' command) by that distribution.
    
    # What you can do is set this up as a loop, where each time you calculate a new photon energy, you then generate a point in the spectrum corresponding to that energy, and then go back to the start, select another random angle/direction, and go again.
    
    # The equations being used below are mostly those defined in R&L, such as the doppler formulae and the formula to calculate the photon energies. Please just be aware that your process may be different, this is just an example of how you might want to implement the equations and calculate the final energy:
    
    # So here is an example of how to start a loop of photons :
    
    spec=[]   # This would be some array that you can save values to
    nphot = int(raw_input("Number of photons 'N_phot': \t")) # This prompt allows you to select how many photons you want
    for i in range(nphot): # Begin loop of photons
      e = 0.01 # units of keV, just as an example of an initial photon energy (notice it's an X-ray)
      e_max = 1e9
      e_min = 0.01
        loge_max = log(e_max) # to get log limits, for binning
            loge_min = log(e_min)


    # Here is just an example of a binning process - it's certainly not directly what you should use, but take a look at this to get a feel for how to set up a logical bin range to dump values to. The author here is generating an electron energy distribution, saving it to an array, and then plotting that distribution. You can also see here which commands you need to set up a plot and save it. Note also that here that this process allows you to choose from a distribution, whereas you will be generating one when you calculate the energies of these photons, which is a bit different :

    ## Generate electron distribution:
    nbinse = 100		# Number of bins for electron energy distribution
    binliste = ([(x+1.0)*(gamma_max-gamma_min)/nbinse for x in range(nbinse)])		# Bins within electron energy domain
    binlistde = [x + gamma_min for x in binliste]	         	# Offset to gamma_minimum
    countse = [int(round(pow(x,2)*sqrt(1.0-1.0/pow(x, 2.0))*exp(- x/(theta_e))/(theta_e*kv(2,1./theta_e))*1e5)) for x in binlistde]	# Counts according to Maxwell-Juttner distribution

    mjdata = []		# Create data list for output electron energy distribution
    for i in range(len(binlistde)):
	x = [binlistde[i],countse[i]]
	mjdata.append(x)

    diste = []
    for i in range(len(binlistde)):		# Electron energy distribution - this would allow one to take a value form the distribution at a later stage - instead you will want to assign a photon energy and counts to a bin
	x = [binlistde[i]]*(countse[i])
	diste = x + diste

    np.savetxt("Data_hea_comp/MJ_spec.dat", mjdata)		# Save Maxwell-Juttner distribution data in 'Data_hea_comp' folder - n.b., defined by user

    loglog(binlistde, countse)		# loglog-plot of electron energy distribution
    xlabel('electron energy [gamma]')
    ylabel('counts')
    title('Maxwell-Juttner distribution for relativistic electrons')
    savefig("Data_hea_comp/MJ_spec_plot.pdf")		# Save plot in 'Data_hea_comp' folder - n.b., defined by user
    show()		# Show plot on screen


    # So this could be a way to start your loop. Afterwards you can randomly select an angle for the photon (as shown above), and then use the formulae in R&L to calculate the final photon energy in the lab frame, and you will need to use a binning technique to add up counts at certain energy bins (you decide what sizes these bins should be). Here is some code you could use to calculate a final energy from scattering - please bear in mind you will still need to go through the calculations and work this out for yourself, or the code will make no sense. Basically it shows the angle/energy transformations between frames, and the result is the final photon energy in the lab frame.

    # BEGIN SCATTERING #

    # K = lab frame, K' = electron rest frame

    # Generate random photon incident angle from dist. P(mu)dmu=(1-beta*mu)dmu (K: muth = cos_theta, K': muth_p = cos_theta_p, '_p'='prime')
    rho_min = pow((1.0-beta), 2.0)/(2.0*beta)
    rho_max = pow((1.0+beta), 2.0)/(2.0*beta)
    rho = random()*(rho_max-rho_min)+rho_min
    muth = (1.0-sqrt(2.0*beta*rho))/beta
    muth_p = (muth-beta)/(1.0-beta*muth)		# Aberration of light: eq. 4.8b, for muth_p

    e1_p = e*gamma*(1.0-beta*muth)		#Final photon energy in electron rest frame:  eq. 7.7 & Thomson approx: e1_p ~ e_ p (elastic scattering)

    # Rejection method for scattering angle mua_p in K' - this just assumes some random direction for the scattering and confines the photon to that region:
    a = 2.0*random()-1.
    b = 2.0*random()
    if b >= 1.0+pow(a,2):
    a = -1.0+2.0*random()
    b = 2.0*random()
    mua_p = a

    muphi_p = cos(random()*2.0*pi)		# Uniform distribution for muphi_p
    sin_theta_p = sqrt(1.0-pow(muth_p,2))		# outgoing theta prime
    sin_alpha_p = sqrt(1.0-pow(mua_p,2))		# outgoing alpha prime
    muth1_p = muth_p*mua_p-sin_alpha_p*sin_theta_p*muphi_p		# scattering angle in electron rest frame: muth1_p = cos_theta1_p

                e1 = e1_p*gamma*(1.0+beta*muth1_p)		# Final photon energy in lab frame



