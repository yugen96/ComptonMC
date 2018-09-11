// Standard libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
// Custom libraries
#include "utils.hh" // Some utility functions
#include "photon.hh" // Photon class
#include "spectrum.hh"
#include "constantsCGS.hh" // Physical and astronomical constants


#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

// Seed for random number generation
std::random_device rd ;
std::mt19937 gen(rd()) ; 
std::uniform_real_distribution<> dis(0,1) ;


///////////// Constructor /////////////
// Specifies the number of photons that constitute the spectrum
// and initializes the energy bin vector
Spectrum::Spectrum(int nr_Y, std::vector<double> bins) : 
        N(nr_Y), Ebins(bins.size()), counts(bins.size(),0), 
        input_counts(bins.size(),0), gammas(nr_Y,0) { 
   
    // Initialize the energy bin array
    int i ;
    for(i=0 ; i<bins.size() ; i++) { Ebins[i] = bins[i] ; }
}
        
        
///////////// Custom functions /////////////
// Generates the cosine of a random angle distributed according to 
// f(mu) = (1 - beta*mu)/2 via the acceptance-rejection method, where mu = cos(theta)
double Spectrum::randMu(double beta) {
    // Draw a random probability from the CDF
    double xi, xi2, fmu ;
    bool found = false ;
    while (!found) {
        xi = 2*dis(gen) - 1 ;
        xi2 = dis(gen) ;
        fmu = 0.5 * (1 - beta*xi) ;
        if (xi2 <= fmu) { found = true ;}
    }
    //std::cout << "In randMu: xi=mu=" << xi << std::endl ; TODO REMOVE
    return xi ;
}




// Calculates a single Compton up-scatter 
// N.B.: this function assumes all variables are given in the scattering electron rest frame
//       (imagine primes behind each variable in the following code if you will)
void Spectrum::comptscat(Photon &phot, double beta) {  
    // Generate a random scattering angle
    double mu_theta1, temp, mu_THETA1, sin_theta, sin_theta1, mu_phi1=dis(gen)*2*M_PI ;
    while(temp >= 1+pow(mu_theta1,2)) {
        mu_theta1 = 2*dis(gen) - 1 ;
        temp = 2*dis(gen) ;
    }
    sin_theta = 1 - pow(phot.getCosTheta(),2) ;
    sin_theta1 = 1 - pow(mu_theta1,2) ;
    mu_THETA1 = phot.getCosTheta()*mu_theta1 - sin_theta*sin_theta1*mu_phi1 ;

    // Compute the upscattered energy
    phot.setMu_th(mu_THETA1) ; 
} 


// The Bose-Einstein distribution
double Spectrum::BE_distr(double f, double T) { 
    return (1.0 / (exp(const_cgs::h*f / (const_cgs::kB*T)) - 1.0)) ; 
}
   
     
// Gives the black-body intensity for temperature T and frequency f
double Spectrum::BB(double f, double T) {
    double dos = 2*const_cgs::h*pow(f,3) / pow(const_cgs::c,2) ;
    return (dos * BE_distr(f,T)) ;
}


// Samples a Black Body function according to the method specified in Carter and Cashwell (1975)
double Spectrum::sampleBB(double T) {
    
    // Find the minimum value l for which sum_1^l 1/j^4 >= xi0 * pi^4 / 90
    int j=1 ;
    bool found = false ;
    double L, temp=0, xi=dis(gen) ;
    while (!found) {
        temp += pow(j,-4) ;
        if (temp >= xi*pow(M_PI,4)/90) { 
            L = j ;
            found = true ; 
        }
        else { j++ ; }
    }
    
    // Compute and return the frequency via x = h*nu/(k_B*T) = -ln(xi1*xi2*xi3*xi4)/L
    int i ;
    temp = 0;
    xi = dis(gen) ;
    for (i=0 ; i<4 ; i++) { 
        xi = dis(gen) ;
        temp += log(xi) ;
    }
    return -(const_cgs::kB * T / L) * temp / const_cgs::h ;
}


// Samples a frequency from the single-particle synchrotron emission 
// spectrum (see R&L (6.31c)) via the acceptance-rejection method 
// Fmax is the maximum value of F(x)
double Spectrum::sampleSynch(double Fmax, double gamma, double B) {
    // Draw a random probability from the CDF
    double xi, xi2, Fx, w_c;
    bool found = false ;
    while (!found) {
        xi = 6 * dis(gen)  ; // See R&L fig.6.6 --> F(x) ~ 0 for x > 4
        xi2 = Fmax * dis(gen) ;
        Fx = xi * utils::int_bessel_k(xi, 6, 5/3.0, 100) ; // Integrate out to 10*w_c
        if (xi2 <= Fx) { found = true ;}
    }
    w_c = ( 3/2.0 * pow(gamma,2) * const_cgs::e * B / 
                            (const_cgs::m_e * const_cgs::c) ) ;
    return w_c * xi / (2*M_PI) ;
}


// Samples a photon energy from the distribution p(E) = 1/E (to be used for the LPA)
// See appendix A of Joerns notes
double Spectrum::sample_finv(double fmin, double fmax) {
    double xi = dis(gen), lgthm = log(fmax/fmin) ;
    return fmin * exp(xi * lgthm) ;
}


// Generates a random electron Lorentz-boost factor from a power-law
// spectrum of electrons with spectral index p via the inversion method
double Spectrum::gen_gamma(double p) {
    double cdf = dis(gen) ; // CDF random value
    return pow(1-cdf,1/(1-p)) ;
}

    
// Finds to which bin a photon belongs and adds its escaped fraction to the bin
bool Spectrum::addY_to_bin(Photon& phot, std::vector<double>& countsvec,
                           double opt_depth, bool LPA) {
    // Find the energy bin to which the photon belongs 
    int i ; 
    bool bin_found = false ; 
    for (i=1 ; i<countsvec.size() ; i++){
        // Check boundaries
        if ((phot.getE() < Ebins[0]) || (phot.getE() > Ebins[Ebins.size()-1])) {
            return false ;
        }
        // Intermediate bins
        else if ((phot.getE() >= Ebins[i-1]) && (phot.getE() < Ebins[i])) {
            bin_found = true ; 
            if (LPA) { countsvec[i] += phot.getw() * exp(-opt_depth) ; }
            else { 
                countsvec[i] += phot.getw() ; 
            }
            return true ;
        }
    }
}


// Simulates the synchrotron part of the spectrum
void Spectrum::sim_synch(double GAMMA, double mu_alpha, double B) {
    
    // Generate N synchrotron photons and put them in their respective energy bins
    bool binfound ;
    double gamma, beta, phot_f, phot_E, phot_E__JRF, BETA ;
    for (int i=0; i<N; i++) {   
        std::cout << i << std::endl ;
        //　Pick jet rest frame electron
        gamma = gen_gamma() ;
        gammas[i] = gamma ;
        beta = sqrt(1-pow(gamma,-2)) ;
        
        // Sample a random photon frequency
        phot_f = sampleSynch(0.915,gamma,B) ; // Hz
        phot_E__JRF = const_cgs::h * phot_f ; // erg

        // Boost back to the observer frame
        BETA = sqrt(1-1/pow(GAMMA,2)) ;
        phot_E = phot_E__JRF * GAMMA * (1 + BETA*mu_alpha) ;
        
        // Create photon object with random angle compared to potential scattering partner
        Photon photon(phot_E, randMu(beta)) ; 
        photon.assign_weight(1,1) ; 
        
        // Add the photon to its corresponding energy bin
        binfound = addY_to_bin(photon,input_counts) ; 
        // If the bin cannot be found, generate a new photon
        if (!binfound) { i-- ; }  
    }
}


// Simulates the spectrum for a jet emission region of size R at temperature T 
// with optical depth tau which has a bulk lorentz factor GAMMA and makes an angle
// alpha with the line of sight such that mu_alpha=cos(alpha).
// LPA can be specified to turn the large particle approach on or off
void Spectrum::sim_comptspect(double T, double tau, double GAMMA,
                              double mu_alpha, bool LPA) {  
    
    // Generate a loop of photons to be scattered
    int i ; // iterator
    bool binfound, scatter; 
    double phot_f, phot_E, phot_E__JRF, n_Y, x, P, gamma, beta, BETA ;
    for (i=0; i<N; i++) {   
        std::cout << i << std::endl ;
        
        // Sample a random photon frequency
        if (LPA) { phot_f = sample_finv(1e13, 1e25) ; } 
        else { phot_f = sampleBB(T) ; }
        phot_E = const_cgs::h * phot_f ; // The photon energy // erg
        
        // Create an input photon object
        Photon photon(phot_E,mu_alpha) ;
        binfound = addY_to_bin(photon,input_counts,tau,false) ; // Add the input photon
        if (!binfound) {
            i-- ;
            continue ;
        }
        
        // Boost photon energy to the jet rest frame (JRF)
        BETA = sqrt(1-1/pow(GAMMA,2)) ;
        phot_E__JRF = phot_E * GAMMA * (1 - BETA*photon.getCosTheta()) ;
        //　Pick jet rest frame electron for potential scattering
        gamma = gen_gamma() ;
        gammas[i] = gamma ;
        beta = sqrt(1-pow(gamma,-2)) ;
        // Create photon object with random angle compared to potential scattering partner
        photon.setE(phot_E__JRF) ;
        photon.setMu_th(randMu(beta)) ; 
        
        // Assign statistical weight to the photon
        if (LPA) { 
            // The input photon number density
            n_Y = 4*M_PI * BB(phot_f,370) / (const_cgs::h * phot_E) ; // # cm-2 s-1 erg-1
            photon.assign_weight(n_Y, 1/phot_E) ; 
        }
        else { photon.assign_weight(1,1) ; }
        
        // Assign escaped photon fraction to correct energy bin
        scatter = true ;
        while (scatter==true) {
            // In case of the large particle appraoch add photon weight to correct
            // photon energy bin at the start of each scattering
            if (LPA) { 
                binfound = addY_to_bin(photon,counts,tau,LPA) ; 
                // Generate new photon if the photon didn't fit into any of the 
                // prespecified energy bins
                if (!binfound) { 
                    i-- ; 
                    break ;
                } 
            }
            
            // Decide whether the photon scatters or not
            x = dis(gen) ; 
            if (photon.getE() > const_cgs::eV) { 
                P = 0.5*pow(const_cgs::eV/photon.getE(),3) ; 
            }
            else { P = 0.5 ; } 
        
            // Compton scatter the photon if scatter==TRUE
            if (x>=P) { scatter = false ; } 
            else {                
                photon.Lorentzboost(beta) ; // Boost to the electron rest frame
                comptscat(photon, beta) ; // Upscatter the photon in the rest frame
                photon.Lorentzboost(-beta) ; // Boost back to the lab frame
            }
        }
        
        // Boost back out of the jet rest frame
        phot_E = photon.getE() * GAMMA * (1 + BETA*mu_alpha) ;
        photon.setE(phot_E) ;

        // In case one doesn't use the large particle approach, just add 1 to the 
        // correct photon energy bin
        if (!LPA) { 
            binfound = addY_to_bin(photon,counts,tau,LPA) ; 
            if (!binfound) { i-- ; }
        }
    }
}




// Writes the spectrum away to a text file
void Spectrum::spectrumtxt(std::string inputY_f, std::string comptY_f, 
                           std::string e_f) const {
    // Open text file
    std::ofstream f_inputY(inputY_f), f_comptY(comptY_f), f_e(e_f) ;
    f_inputY << "#Input photon Spectrum\n-------------------" << std::endl ;
    f_inputY << "#Bin\t \t#counts\n----        -------" << std::endl ;
    // Write input photon distribution data to file
    int i;
    for (i=0 ; i<Ebins.size() ; i++) { 
        f_inputY << Ebins[i] << "\t,\t" << input_counts[i]/N << std::endl ; 
    }
    f_inputY.close() ;

    // Open text file
    f_comptY << "#Compton scattered photon Spectrum\n-------------------" << std::endl ;
    f_comptY << "#Bin\t \t#counts\n----        -------" << std::endl ;
    // Write photon distribution data to file
    for (i=0 ; i<Ebins.size() ; i++) { 
        f_comptY << Ebins[i] << "\t,\t" << counts[i]/N << std::endl ; 
    }
    f_comptY.close() ;
    
    // Write electron distribution data to file
    f_e << "#Electron Lorentz Factors\n-------------------" << std::endl ;
    f_e << "#Electron Number\t \t#Lorentz factor" << std::endl ;
    f_e << "----------------\t \t---------------" << std::endl ;
    // Write data to file
    for (i=0 ; i<N ; i++) { 
        f_e << i << "\t\t\t\t\t,\t" << gammas[i] << std::endl ; 
    }
    // Close file
    f_e.close() ;
}


// Access functions
std::vector<double> Spectrum::getEbins() const { return Ebins ; } 
std::vector<double> Spectrum::getcounts() const { return counts ; } 






