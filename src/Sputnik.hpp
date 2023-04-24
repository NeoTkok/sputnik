#ifndef SPUTNIK_HPP
#define SPUTNIK_HPP

#include "Orbit.hpp"

class sputnik {
private:
    orbit m_orb;
    double g;
public:
    sputnik(const orbit A, const double G) : m_orb(A), g(G) {};

    double get_h(){
        return - g * (1-m_orb.e * m_orb.e) / m_orb.getP();
    }

// радиальная скорость
    double v_r(double alpha){
        return sqrt(g/m_orb.getP()) * m_orb.e * sin(alpha);
    }

// трансверсальная скорость
    double v_tau(double alpha){
        return sqrt(g / m_orb.getP()) * (1 + m_orb.e * cos(alpha));
    }

//***********
// скорость в перицентре
    double v_p(){
        return (1 + m_orb.e) * sqrt(g / m_orb.getP());
    }

//***************
// новый фокальный параметр при изменении скорости в перицентре на dv
    double new_p(double dv){
        double c = (v_p() + dv) * m_orb.getP() / (1 + m_orb.e);
        return c * c / g;
    }


//****************
// новый  при изменении скорости в перицентре на dv
    double new_e(double dv){
        return (v_p() + dv) * (v_p() + dv) * m_orb.getP() / (1 + m_orb.e) / g - 1;
    }

    double get_v(double phi){
        return sqrt(get_h() + 2. * g / m_orb.getR(phi)) ;
    }


    double manevr(orbit A){
        sputnik B(A, g);

        std::vector<double> Z = null_f(m_orb, A);
        double ugol = exact_null_f(m_orb, A, Z[1], Z[2], 1e-8);
        
        double b = cos_beta(m_orb, A, ugol);
        double v0 = get_v(ugol);
        double v = B.get_v(ugol);

        return sqrt(v*v + v0*v0 - 2*v*v0*b);
    
    }

};


#endif 