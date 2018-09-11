// Standard libraries
#include <iostream>
#include <ctime>
#include <cmath> 
#include <vector> 
// Custom libraries
#include "constantsCGS.hh" // Physical and astronomical constants
#include "photon.hh"



///////////// Constructor /////////////
// Creates a photon with specified energy and random direction vector
Photon::Photon(double energy, double costheta, double W) : 
        E(energy), mu_th(costheta), w(W) {}       


///////////// Custom functions /////////////
// Determines the relativistic angle abberation when transforming from lab to rest frame
double Photon::transf_mu(double mu, double beta) { return (mu-beta) / (1-beta*mu) ; }

// Conducts a Lorentz translation along the z-direction
void Photon::Lorentzboost(double beta) {
    // Determine new energy
    double gamma = 1 / sqrt(1-pow(beta,2)) ;
    double new_E = gamma*E*(1 - beta*mu_th) ;
    E = new_E ;
    // Determine new direction
    mu_th = transf_mu(mu_th,beta) ; 
}


///////////// Access functions ///////////      
// Assign weight to the photon based on the input photon number density 
// and energy probability distribution p(E) ~ E^(-1)
void Photon::assign_weight(double N_ph, double p_E) { w = N_ph*p_E ; }
// Assign energy
void Photon::setE(double energy) { E = energy ; }
// Assign direction
void Photon::setMu_th(double costheta) { mu_th = costheta ; }

///////////// Read functions /////////////
double Photon::getE() const { return E ; }
double Photon::getw() const { return w ; }
double Photon::getCosTheta() const { return mu_th ; }

///////////// Print functions /////////////
void Photon::printE() const { std::cout << E << std::endl ; }       









