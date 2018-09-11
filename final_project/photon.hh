#ifndef PHOTON_HH
#define PHOTON_HH

#include "constantsCGS.hh"


class Photon {
    public:
        Photon(double energy=1e2*const_cgs::eV, 
               double costheta=0, double W=1) ;     
        
        double transf_mu(double mu, double beta) ;
        void Lorentzboost(double beta) ;
        void comptscat(double beta) ;

        void setE(double energy) ;
        void assign_weight(double N_ph, double p_E) ;
        void setMu_th(double costheta) ;
        
        double getE() const ;
        double getw() const ;
        double getCosTheta() const ;
        
        void printE() const ;
        
    
    private:
        double E ; // Photon energy
        double w ; // Photon weight
        double mu_th ; // The cosine of the longitudinal angle (mu_th=cos(theta))
        //std::vector<double> a(3); // Cartesian photon propagation direction (x,y,z)
} ;                                     

#endif 
