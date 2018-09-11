#ifndef CONSTANTSCGS_HH
#define CONSTANTSCGS_HH

#include <math.h>

namespace const_cgs {
    // Physical constants
    static const double kB = 1.380658e-16 ; // erg K-1
    static const double c = 2.99792458e10 ; // cm s-1
    static const double c_s = pow(10,7.5) ; // cm s-1
    static const double h = 6.6260755e-27 ; // erg s
    static const double hbar = h / (2*M_PI) ; // erg s
    static const double G = 6.67259e-8 ; // cm3 g-1 s-2
    static const double e = 4.8032068e-10 ; // esu
    static const double m_e = 9.1093897e-28 ; // g
    static const double m_p = 1.6726231e-24 ; // g
    static const double m_n = 1.6749286e-24 ; // g
    static const double amu = 1.6605402e-24 ; // g
    static const double N_A = 6.0221367e23 ; // --
    static const double eV = 1.6021772e-12 ; // erg
    static const double E_e = 0.5e6 * eV ; // erg 
    static const double sigma_sb = 5.67051e-5 ; // erg cm-2 K-4 s-1
    static const double sigma_Th = 6.6524e-25 ; // cm2

    // Astronomical constants
    static const double AU = 1.49597871e13 ; // cm
    static const double pc = 3.08567758e18 ; // cm
    static const double lyr = 9.4605284e17 ; // cm
    static const double Msol = 1.9891e33 ; // g
    static const double Rsol = 6.957e11 ; // cm
    static const double Lsol = 3.839e33 ; // erg s-1
    static const double yr = 365.25*24*3600 ; // s
    static const double H0 = 71e5/(1e6*pc) ; // s-1
}


#endif 
