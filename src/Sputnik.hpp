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
        return - g * (1-m_orb.e * m_orb.e) / m_orb.p;
    }

// радиальная скорость
    double v_r(double alpha){
        return sqrt(g/m_orb.p) * m_orb.e * sin(alpha);
    }

// трансверсальная скорость
    double v_tau(double alpha){
        return sqrt(g / m_orb.p) * (1 + m_orb.e * cos(alpha));
    }

//***********
// скорость в перицентре
    double v_p(){
        return (1 + m_orb.e) * sqrt(g / m_orb.p);
    }

//***************
// новый фокальный параметр при изменении скорости в перицентре на dv
    double new_p(double dv){
        double c = (v_p() + dv) * m_orb.p / (1 + m_orb.e);
        return c * c / g;
    }


//****************
// новый  при изменении скорости в перицентре на dv
    double new_e(double dv){
        return (v_p() + dv) * (v_p() + dv) * m_orb.p / (1 + m_orb.e) / g - 1;
    }

// получение скорости при заданном угле
    double get_v(double phi){
        return sqrt(get_h() + 2. * g / m_orb.getR(phi)) ;

// получение характерной скорости при заданном угле(в котором происходит пересечение орбит)
    }
    double delta_V(orbit A, double phi){
        sputnik B(A,g);
        double b = cos_beta(m_orb, A, phi); // угол межнду начальной скоростью и целевой
        double v0 = get_v(phi); // как раз таки начальная скорость
        double v = B.get_v(phi); // скорость, которую будет иметь КА на целевой орбите

        return sqrt(v*v + v0*v0 - 2*v*v0*b); // т. Косинусов
    }
//****** 
// получение скоростей при ОЧЕНЬ малых углах(где касательные)
    double delta_Vr(orbit A, double phi){
        sputnik B(A,g);
        double v0 = get_v(phi); 
        double v = B.get_v(phi); 
        return v - v0;
    }


//
    double manevr(orbit A, double eps){ 
        sputnik B(A, g); 
        std::array<double,3> a = peresech(m_orb, A); // направляющий вектор пересечения орбит

        if(a[0] == 0 && a[1] == 0 && a[2] == 0) // т.е. лежит в плоскости
        {
            std::cout<< "плоскость целевого эллипса совпала с настоящей" << std::endl;

            std::vector<double> x = intersection(m_orb, A, eps); // экстремальные точке
            std::vector<double> V; // вектор скоростей
            std::vector<double> PHI; // вектор углов из которых мы видим точку пересечения
            std::vector<double> UGOL; // угол между скоростями(начальной и целевой)

            if (x.size() == 1){ // если экстремальные точки 
                V.push_back(delta_Vr(A, x[0]));
                PHI.push_back(x[0]);
                UGOL.push_back(0);
            }
            if (x.size() == 2){  // добавление двух элементов
                V.push_back(delta_V(A, x[0]));
                PHI.push_back(x[0]);
                UGOL.push_back(cos_beta(m_orb, A, x[0]));
                V.push_back(delta_V(A, x[1]));
                PHI.push_back(x[1]);
                UGOL.push_back(cos_beta(m_orb, A, x[1]));
            }

            if (x.size() == 0){
            }


        for (int i = 0; i < V.size(); ++i)
            std::cout << UGOL[i] << " - " << V[i] << " - " << PHI[i] << std::endl;
    // таким образом мы получили 3 вектора из которых можно найти оптимальную скорость
    // если эти векторы пустые, то пересечений вообще нет => изменять скорость в перицентре
        }

        
        return -1;
    }

};





#endif 