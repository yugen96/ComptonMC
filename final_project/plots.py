import numpy as np 
import matplotlib.pyplot as plt 


# Constants
h = 6.6260755e-27 # erg s
m_e = 9.1093897e-28 # g
c = 2.99792458e10 # cm s-1
eV = 1.6021772e-12 # erg
E_e = 0.5e6 * eV # erg 
sigma_Th = 6.6524e-25 # cm2
pc = 3.08567758e18


# Source parameters
#----------- jet parameters -----------#
mu_alpha = np.cos(6*np.pi/180) # --
GAMMA =15 # -- 
BETA = np.sqrt(1-GAMMA**(-2)) # --
delta = 1/(GAMMA*(1-BETA*mu_alpha)) # --
#----------- distance parameters -----------#
H0 = 71e5/(1e6*pc) # cm s-1 cm-1
z = 0.859 # -- 
d = c*z/H0 # cm
R = 5e16 # cm
#----------- electron distribution parameters -----------#
gamma2_e__avg = 769.0 # --
n_e = 1 # cm-3
N_e = n_e * (4/3.0) * np.pi * R**3 # --
#----------- radiation density parameters -----------#
U_phot = 5e3 # erg cm-3


# Factor to convert from the photon histogram to a SED
#prefact = n_e * 4/3 * sigma_Th * c * gamma2_e__avg * delta**4 * GAMMA**2 * U_phot * R**3 * 1e6


# Read in SED data
freq, nu_Fnu, sigma_nu_Fnu = np.loadtxt("3C454.3_SEDnuFnu.txt",comments=("#"), usecols=(0,2,3)).T


# Read in the photon and electron energy distribution data
Ebins, inputPhot_counts = np.loadtxt("spectrum_inputY.txt", comments=("#","--"), delimiter=",").T
Ebins, comptPhot_counts = np.loadtxt("spectrum_comptY.txt", comments=("#","--"), delimiter=",").T
# Read in the electron distribution data
e_gammas = np.loadtxt("spectrum_e.txt", comments=("#","--"), delimiter=",", 
                      usecols=1, skiprows=int(5e5)).T 
                      # Only take 5e5 out of 1e6 electrons to speed up the loading


# The electron distribution
fig = plt.figure()
plt.hist(e_gammas, np.logspace(0,4,1000), log=True, normed=1)
plt.xscale('log')
plt.xlabel(r"$\gamma_e$ [--]", fontsize=16)
plt.ylabel(r"normalized counts [--]", fontsize=16)
plt.tight_layout()
plt.savefig("3C454_3_electrondist.png")
#plt.show()
plt.close()

# The initial photon distribution
fig = plt.figure()
ax = plt.axes()
ax.bar((Ebins/h)[1::], inputPhot_counts[1::]*(Ebins[1::]), np.diff(Ebins/h))
ax.set_xscale("log", nonposx='clip'), ax.set_yscale("log", nonposy='clip')
ax.set_xlim(xmax=1e15)
ax.set_xlabel(r"$\nu$ [Hz]", fontsize=16)
ax.set_ylabel(r"$\nu F_{\nu}$ [erg cm-2 s-1]", fontsize=16)
plt.tight_layout()
plt.savefig("3C454_3_inputSED.png")
#plt.show()
plt.close()

# The upscattered photon distribution
fig = plt.figure()
ax = plt.axes()
ax.bar((Ebins/h)[1::], comptPhot_counts[1::]*(Ebins[1::]), np.diff(Ebins/h))
ax.errorbar(10**freq, 10**nu_Fnu, fmt='o', markersize=1, capsize=0.05, color='k')
ax.set_xscale("log", nonposx='clip'), ax.set_yscale("log", nonposy='clip')
#plt.errorbar(10**freq, 10**nu_Fnu, yerr=10**sigma_nu_Fnu, color='0.75', fmt='o', markersize=1)
ax.set_xlabel(r"$\nu$ [Hz]", fontsize=16)
ax.set_ylabel(r"$\nu F_{\nu}$ [erg cm-2 s-1]", fontsize=16)
plt.tight_layout()
plt.savefig("3C454_3_comptSED.png")
#plt.show()
plt.close()


'''
inputPhot_countsLPA = np.loadtxt("spectrum_inputY__LPA.txt", comments=("#","--"), delimiter=",", usecols=1).T
comptPhot_countsLPA = np.loadtxt("spectrum_comptY__LPA.txt", comments=("#","--"), delimiter=",", usecols=1).T
phot_countsSYNCH = np.loadtxt("spectrum_inputY__SYNCH.txt", comments=("#","--"), delimiter=",", usecols=1).T


fig = plt.figure()
ax = plt.axes()
ax.bar((Ebins/h)[1::], spectrum_inputY__LPA[1::]*(Ebins[1::]**2), np.diff(Ebins/h))
ax.errorbar(10**freq, 10**nu_Fnu, fmt='o', markersize=1, capsize=0.05, color='k')
ax.set_xscale("log", nonposx='clip'), ax.set_yscale("log", nonposy='clip')
#plt.errorbar(10**freq, 10**nu_Fnu, yerr=10**sigma_nu_Fnu, color='0.75', fmt='o', markersize=1)
ax.set_xlabel(r"$\nu$ [Hz]", fontsize=16)
ax.set_ylabel(r"$\nu F_{\nu}$ [erg cm-2 s-1]", fontsize=16)
ax.set_ylim(ymin=1e-20)
plt.tight_layout()
plt.savefig("3C454.3_SED__LPA.png")
#plt.show()
plt.close()

fig = plt.figure()
ax = plt.axes()
ax.bar((Ebins/h)[1::], spectrum_comptY__LPA[1::]*(Ebins[1::]**2), np.diff(Ebins/h))
ax.errorbar(10**freq, 10**nu_Fnu, fmt='o', markersize=1, capsize=0.05, color='k')
ax.set_xscale("log", nonposx='clip'), ax.set_yscale("log", nonposy='clip')
#plt.errorbar(10**freq, 10**nu_Fnu, yerr=10**sigma_nu_Fnu, color='0.75', fmt='o', markersize=1)
ax.set_xlabel(r"$\nu$ [Hz]", fontsize=16)
ax.set_ylabel(r"$\nu F_{\nu}$ [erg cm-2 s-1]", fontsize=16)
ax.set_ylim(ymin=1e-20)
plt.tight_layout()
plt.savefig("3C454.3_SED__LPA.png")
#plt.show()
plt.close()


fig = plt.figure()
ax = plt.axes()
ax.bar((Ebins/h)[1::], spectrum_inputY__SYNCH[1::]*(Ebins[1::]**2), np.diff(Ebins/h))
ax.errorbar(10**freq, 10**nu_Fnu, fmt='o', markersize=1, capsize=0.05, color='k')
ax.set_xscale("log", nonposx='clip'), ax.set_yscale("log", nonposy='clip')
#plt.errorbar(10**freq, 10**nu_Fnu, yerr=10**sigma_nu_Fnu, color='0.75', fmt='o', markersize=1)
ax.set_xlabel(r"$\nu$ [Hz]", fontsize=16)
ax.set_ylabel(r"$\nu F_{\nu}$ [erg cm-2 s-1]", fontsize=16)
ax.set_ylim(ymin=1e-20)
plt.tight_layout()
plt.savefig("3C454.3_SED__SYNCH.png")
plt.show()
plt.close()
'''



