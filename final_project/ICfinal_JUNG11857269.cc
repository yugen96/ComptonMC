// Standard libraries
#include <iostream> // in- and output
#include <cmath> // Basic arithmetic
#include <vector> // Standard vector library
#include <random> // For random number generation
#include <string> 

// Some libraries and classes
#include "utils.hh"
#include "photon.hh" // Photon class
#include "spectrum.hh" // Spectrum class
#include "constantsCGS.hh" // Physical and astronomical constants



int main() {	

    // Source parameters
    //------- temperature and magnetic field -------//
    double T = 340 ; // Source temperature // K
    double B = 1 ; // Magnetic field strength inside jet // T
    //------- Accretion info -------//
    double M_BH = 2.3e9 * const_cgs::Msol ; // Central BH mass // g 
    double L_edd = 1e47 ; // Eddington Luminosity // erg s-1
    double Mdot_edd = L_edd / (0.3 * pow(const_cgs::c,2)) ; // Eddington accretion rate assuming
                                                            // efficiency of 30% // g s-1
    //------- distances -------//
    double R_g = 2*M_BH*const_cgs::G / pow(const_cgs::c,2) ; // The Schwarzschild radius // cm
    double R_em = 5e16 ; // The emission region radius // cm
    double z = 0.859 ; // The redshift // --
    double d = z * const_cgs::c / const_cgs::H0 ; // Distance to the source // cm
    //------- jet parameters -------//
    double GAMMA = 15 ; // Jet Bulk factor // --
    double alpha = 6/180*M_PI ; // Jet inclination angle // rad
    //------- particle distribution -------//
    // An estimate for the particle number density in the accretion region // cm-3
    double n = Mdot_edd / (4*M_PI*pow(R_g,2) * const_cgs::c_s * const_cgs::m_p) ; //TODO REMOVE
    double n_e = 1 ; // Assume there are of the order unity electrons cm-3 //TODO DEFEND!!
    double tau = n_e * const_cgs::sigma_Th * R_em ; // The associated optical depth // --
    

    // Define the energy bins
    std::vector<double> E_bins = utils::logspace(13,25,100) ;
    utils::vecmult(const_cgs::h, E_bins) ;
    // Create the compton scattered photon spectrum
    Spectrum spectr(1e6, E_bins) ;
    spectr.sim_comptspect(T, 0.1, GAMMA, cos(alpha), false) ;
    spectr.spectrumtxt("spectrum_inputY2.txt", "spectrum_comptY2.txt", "spectrum_e2.txt") ;
    
    
    /* /TODO Large Particle Approximation
    // Define the energy bins
    std::vector<double> E_bins = utils::logspace(13,25,100) ;
    utils::vecmult(const_cgs::h, E_bins) ;
    // Create the compton scattered photon spectrum
    Spectrum spectr(1e6, E_bins) ;
    spectr.sim_comptspect(T, 0.1, GAMMA, cos(alpha), true) ;
    spectr.spectrumtxt("spectrum_inputY__LPA.txt", "spectrum_comptY__LPA.txt", 
                       "spectrum_e__LPA.txt") ;
    */
    
    
    /* /TODO Synchrotron part
    // Create the compton scattered photon spectrum
    std::vector<double> synch_bins = utils::logspace(6,20,100) ;
    utils::vecmult(const_cgs::h, synch_bins) ;
    // Create the synchrotron part
    Spectrum synch_spectr(1e5, synch_bins) ;
    synch_spectr.sim_synch() ;
    synch_spectr.spectrumtxt("spectrum_inputY__SYNCH.txt", "spectrum_comptY__SYNCH.txt", 
                             "spectrum_e__SYNCH.txt") ;
    */
    
    return 0 ;
}



