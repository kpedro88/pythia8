
#ifndef SUEPSHOWER_h
#define SUEPSHOWER_h

#include <iostream>
#include <math.h>
#include <boost/math/tools/roots.hpp>
#include <boost/bind.hpp>
#include <vector>

class tolerance {
public:
	tolerance(double eps) :
		_eps(eps) {
	}
	bool operator()(double a, double b) {
		return (fabs(b - a) <= _eps);
	}
private:
	double _eps;
};



class Suep_shower 
{ 
    public: 
  
    // Data Members 
    double m;
    double Temp;
    double Etot; 
    
    // Constructor
    Suep_shower(double mass, double temperature, double energy);
    
    // methods
    double f(double p);
    double fp(double p);
    double test_fun(double p);
    std::vector<double> generate_fourvector(); 
    std::vector< std::vector <double> > generate_shower();
    double reballance_func(double a, std::vector< std::vector <double> > event);
        
    private:
    // private variables
    double A;
    double p_m;
    std::pair <double,double> p_plusminus;
    double p_plus, p_minus;
    double lambda_plus,lambda_minus;
    double q_plus, q_minus, q_m;
            
}; 

 

#endif





