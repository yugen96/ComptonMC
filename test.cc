


// Standard libraries
#include <iostream>
#include <ctime>
#include <vector>
// Custom libraries
#include "utils.hh" // Some utility functions
#include "photon.hh" // Photon class
#include "constantsCGS.hh" // Physical and astronomical constants
#include "boost/math/interpolators/barycentric_rational.hpp"





int main() {
    
    std::vector<int> intvec(10) ;
    std::vector<int>::iterator iter = intvec.begin() ;
    int i = 10;
    while(iter!=intvec.end()) {
        *iter = i ;
        std::cout << *iter << std::endl ;
        iter ++ ;
        i-=2 ;
    }
    

    return 0 ;
}
