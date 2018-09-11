// Imports
#include <iostream>
#include <algorithm> // std::is_sorted
#include <vector>
#include <cmath>
#include "utils.hh"
//#include "boost/math/special_functions/bessel.hpp"




///////////////// STD::VECTOR RELATED FUNCTIONS /////////////////
// Prints the contents of a vector on the screen
void utils::printvec(const std::vector<double> &v) {
    std::vector<double>::const_iterator iter = v.begin() ;
    while(iter!=v.end()) {
        std::cout << *iter << std::endl ;
        iter++ ;
    }
}

// Multiplies a vector and a scalar
void utils::vecmult(double scalar, std::vector<double> &v) {
    std::vector<double>::iterator iter=v.begin() ;
    while(iter!=v.end()) {
        *iter *= scalar ;
        iter++ ;
    }
}

// Returns a linearly spaced array between begin and end with specified
// stepsize and number of elements
std::vector<double> utils::linspace(double begin, double end, int num) {
    // Define and initialize returnvalue vector
    std::vector<double> retval(num) ;
    retval[0] = begin ;
    
    // Fill the array
    double stepsize = (end - begin) / (num-1) ;
    std::vector<double>::iterator iter = retval.begin()+1 ;
    while(iter!=retval.end()) { 
        *iter = *(iter-1) + stepsize ; 
        iter++ ;
    } 
    return retval ;        
}

// Returns a logarithmically spaced array starting at base^start and stopping at base^stop
// with a total number of steps given by num
std::vector<double> utils::logspace(double begin, double end, int num, double base) {
    // Define and initialize returnvalue vector
    std::vector<double> retval(num) ;
    retval[0] = pow(base,begin) ;
    
    // Fill the array
    double prevpow=begin, currpow ;
    double stepsize = (end - begin) / (num-1) ;
    std::vector<double>::iterator iter = retval.begin()+1 ;
    while(iter!=retval.end()) { 
        currpow = prevpow+stepsize ;
        *iter = pow(base,currpow) ;
        prevpow=currpow ;
        iter++ ;
    } 
    return retval ;        
}

// Zips two vectors
void utils::zip(const std::vector<double> &v1, 
         const std::vector<double> &v2,
         std::vector<std::pair<double,double>> &zipped) {
    for (int i=0 ; i<v1.size() ; i++) {
        zipped.push_back(std::make_pair(v1[i],v2[i])) ;
    }
}    

// Unzips a zipped vector
void utils::unzip(const std::vector<std::pair<double,double>> &zipped,
           std::vector<double> &v1, 
           std::vector<double> &v2) {
    for (int i=0 ; i<v1.size() ; i++) {
        v1[i] = zipped[i].first ;
        v2[i] = zipped[i].second ;
    }
}

// Allows for sorting a vector of pairs by each element's second entry
bool utils::sort_by_second::operator()(std::pair<double,double> pair1, 
                                       std::pair<double,double> pair2) { 
        return (pair1.second < pair2.second) ;
} 

// Returns interpolated value at x from parallel vectors (xData , yData)
// Input parameters: 
//      - x:            the point of interpolation
//      - xData:        the x-data from which to interpolate
//      - yData:        the y-data from which to interpolate
//      - extrapolate:  allows for doing extrapolation at the endpoints
//      - sort:         allows for sorting the x-data if not done already
double utils::linpol(double x, std::vector<double> &xData, std::vector<double> yData, 
              bool extrapolate, bool sort) {
    // Checks whether the x-data is sorted to ascending order
    if ( (!std::is_sorted(xData.begin(),xData.end()))  && (sort==true)) {
        std::vector<std::pair<double,double>> zipped ;
        utils::zip(yData,xData,zipped) ;
        std::sort(zipped.begin(), zipped.end(), utils::sort_by_second()) ;
        unzip(zipped,yData,xData) ;
    }        
    else if (sort==true) {
        throw "In utils::linpol: xData is not sorted to ascending order" ; 
    }
    
    // Finds the left end of the interpolation interval
    int i = 0 ; // index for looping
    int size = xData.size() ; // array size
    if (x >= xData[size - 2]) { // Check whether x lies beyond the domain 
        i = size - 2 ;
    }
    else {
        while(x > xData[i+1]) { i++ ; }
    }
    double xL = xData[i], xR = xData[i+1], yL = yData[i], yR = yData[i+1] ;
    
    // If extrapolation is disabled, take the final array values as output
    if (!extrapolate) {
        if ( x < xL ) yR = yL ; // x lies left of the specified domain
        if ( x > xR ) yL = yR ; // x lies right of the specified domain
    }
    
    // Compute the derivative and the interpolated value at the interpolation point
    double slope = (yR - yL) / (xR - xL) ; 
    return (yL + slope * (x-xL)) ;
}


// Compute the Riemann sum over the modified bessel function of second kind from xmin to xmax 
double utils::int_bessel_k(double xmin, double xmax, 
                           double nu, double N) {
    double left = xmin, width = (xmax-xmin)/N ;
    double total = 0 ;
    std::cout << "Warning! utils::int_bessel_k has been commented out temporarily.\n" ;
    std::cout << "Used 1 to compute integral over Bessel function." << std::endl ;
    for (int i=0 ; i<N ; i++) {
        //TODO UNCOMMENT the line below and remove the 1 in case you want to 
        //TODO look this function.
        total += 1 ; //boost::math::cyl_bessel_k(nu, left+width/2.0) * width ; 
        left += width ;
    }
    return total ;
} 



