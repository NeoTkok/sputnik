#ifndef ORBIT_HPP
#define ORBIT_HPP

#include<cmath>

struct orbit{

    const double a; 
    const double e;
    const double OMEGA;
    const double omega;
    const double i;

    orbit(const double a, const double e, const  double omega, const double psi,
          const double i) : a(a), e(e), OMEGA(OMEGA), omega(omega), i(i) {};


    double getRp() const{
        return a * (1 - e);
    }
    

    double getRa() const{
        return a * (1 + e);        
    }


    double getP() const{
        return a * (1 - e * e);
    }

    std::array<double, 3> getPlane() const{
        std::array<double, 3> r;
        r[0] = sin(i) * sin(omega);
        r[1] = -sin(i) * cos(omega);
        r[2] = cos(i);
        return r;
    }

};

#endif 
