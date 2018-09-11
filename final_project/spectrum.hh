#ifndef SPECTRUM_HH
#define SPECTRUM_HH

#include <string>
#include "photon.hh"
#include "constantsCGS.hh"

class Spectrum {
    public: 
        Spectrum(int nr_Y, std::vector<double> bins) ;
        
        double BE_distr(double f, double T=370) ;
        double BB(double f, double T) ;
        double sample_finv(double fmin, double fmax) ;
        double sampleBB(double T) ;
        double sampleSynch(double Fmax, double gamma=4, double B=10) ;
        double gen_gamma(double p=2.36) ;
        bool addY_to_bin(Photon& phot, std::vector<double>& countsvec,
                         double opt_depth=0.1, bool LPA=false) ;
        
        double randMu(double beta) ;
        void comptscat(Photon& phot, double beta) ;
        
        void sim_comptspect(double T=370, double tau=1e-9, double GAMMA=15,
                            double alpha=cos(6/180*M_PI), bool LPA=false) ;
        void sim_synch(double GAMMA=15, double mu_alpha=cos(6/180*M_PI), double B=1) ;
        
        void spectrumtxt(std::string inputY_f, std::string comptY_f, std::string e_f) const ;
        
        std::vector<double> getEbins() const ;
        std::vector<double> getcounts() const ;
        
    private: 
        int N ; 
        std::vector<double> Ebins ;
        std::vector<double> input_counts ;
        std::vector<double> counts ; 
        std::vector<double> gammas ;
        
} ;

#endif 
