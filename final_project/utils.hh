#ifndef UTILS_HH
#define UTILS_HH

namespace utils{
    
    void printvec(const std::vector<double> &v) ;
    
    void vecmult(double scalar, std::vector<double> &v) ;
    
    std::vector<double> linspace(double begin, double end, int num=50) ;
    
    std::vector<double> logspace(double begin, double end, int num=50, double base=10) ;
    
    void zip(const std::vector<double> &v1, 
             const std::vector<double> &v2, 
             std::vector<std::pair<double,double>> &zipped) ; 
    
    void unzip(const std::vector<std::pair<double,double>> &zipped, 
               std::vector<double> &v1, 
               std::vector<double> &v2) ;
    
    struct sort_by_second {
        bool operator()(std::pair<double,double> pair1, std::pair<double,double> pair2) ;
    } ;
    
    double linpol(double x, std::vector<double> &xData, std::vector<double> yData, 
                  bool extrapolate=true, bool sort=true) ;  

    double int_bessel_k(double xmin, double xmax, 
                       double nu=3/5.0, double N=1000.0) ;
    
} 
#endif 
