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
    double new_p_p(double dv){
        double c = (v_p() + dv) * m_orb.p / (1. + m_orb.e);
        return c * c / g;
    }

//****************
// новый e при изменении скорости в перицентре на dv
    double new_e_p(double dv){
        return (v_p() + dv) * (v_p() + dv) * m_orb.p / (1. + m_orb.e) / g - 1.;
    }

// новый фокальный параметр при изменении скорости в перицентре на dv
    double new_p_a(double dv){
        double c = (v_p() + dv) * m_orb.p / (1. - m_orb.e);
        return c * c / g;
    }
//****************
// новый e при изменении скорости в перицентре на dv
    double new_e_a(double dv){
        return 1 - (v_p() + dv) * (v_p() + dv) * m_orb.p / (1. - m_orb.e) / g ;
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

struct Orbit_and_V{
    orbit ORB;
    double SPEED;
    Orbit_and_V(const orbit ORB, const double SPEED) : ORB(ORB), SPEED(SPEED) {}
};

// скорость которую надо поддать или сбавить в перицентре чтобы произошло касание
Orbit_and_V v_min_p(orbit A, double eps){
    Orbit_and_V W(m_orb,0.);
    if(delta_f(m_orb,A,0.) > 0.){ // целевая орбита внутри
        double rp = m_orb.getRp();
        double Vo = sqrt(g / rp);
        orbit new_A = m_orb;
        new_A.e = 0;
        new_A.p = rp;
        int n = intersection(new_A, A, eps).size();
        if (n == 1) // +целевая орбита касается окружности радиуса перицентр    
        {
            W.ORB = new_A;    
            W.SPEED =  Vo - v_p();
        }
        if (n == 0){ // целевая орбита внутри окрцжности разности перецентр 
            double vl = 0.; // левая граница скорости
            double vr = Vo; // правая граница скорости 
            double V = (vl + vr) / 2.; // серединка скорости
            int j = 0;
            new_A.omega += M_PI;
            if (new_A.omega >= 2*M_PI)
                new_A.omega -= 2 * M_PI;
            sputnik G(new_A, g);
            orbit new_B = new_A; // новая орбита после добавления скорости            
            while(j < 25){ // 2^25 ошибка убывает
                ++j;
                // стоит заметить, что скорость добавляем в перицентре                    
                new_B.e = G.new_e_a(V - Vo); // новый эксцентриситет
                new_B.p = G.new_p_a(V - Vo); // новый фокальный параметр
                int k = intersection(new_B, A, eps).size();
                if (k == 1) // новая орбита идеально касается целевой
                    break;
                if (k == 2) // новая орбита пересекла целевую(переборщили)
                    vl = V;
                if (k == 0) // новая орбита ещё внутри                        vr = V;
                    vr = V;    
                V = (vl + vr) / 2.;
            }
            W.ORB = new_B;
           
            W.SPEED = V - v_p();
            return W;
        }
        if (n == 2){  // пересекает окружность    
            double vl = Vo; // левая граница скорости
            double vr = v_p(); // правая граница скорости
            double V = (vl + vr) / 2.; // серединка скорости 
            int j = 0;
            
            orbit new_B = m_orb; // новая орбита после добавления скорости
            while(j < 25){ // 2^25 ошибка убывает
                ++j;
                // стоит заметить, что скорость добавляем в перицентре
                new_B.e = new_e_p(V - v_p()); // новый эксцентриситет
                new_B.p = new_p_p(V - v_p()); // новый фокальный параметр
                int k = intersection(new_B, A, eps).size();
                if (k == 1) // новая орбита идеально касается целевой
                    break;
                if (k == 2) // новая орбита пересекла целевую(переборщили)
                    vl = V;
                if (k == 0) // новая орбита ещё внутри
                    vr = V;
                V = (vl + vr) / 2.;
                }
            W.ORB = new_B;
            W.ORB.omega += M_PI;
            if (W.ORB.omega >= 2*M_PI)
                W.ORB.omega -= 2 * M_PI;
            W.SPEED =  V - v_p();
        }
    }
    else{   // целевая орбита находится снаружи
        orbit new_A = m_orb;
        double vl = v_p();
        double vr = sqrt((1.87)*g*(1+m_orb.e)/m_orb.p); // иначе ломается
            
        new_A.e = new_e_p(vr - v_p()); // новый эксцентриситет
        new_A.p = new_p_p(vr - v_p()); // новый фокальный параметр            
        int k = intersection(new_A, A, eps).size();
            
        if (k == 0){
            std::cout << "нельзя добавить в перцентр скорость" << std::endl;
            return W;
        }
        double V = (vl + vr) / 2.; // серединка скорости 
        int j = 0;
        while(true){ // 2^15 ошибка убывает
            ++j;
            new_A.e = new_e_p(V - v_p()); // новый эксцентриситет
            new_A.p = new_p_p(V - v_p()); // новый фокальный параметр
            
            int k = intersection(new_A, A, eps).size();
            if (k == 1) // новая орбита идеально касается целевой
                break;
            if (k == 2) // новая орбита пересекла целевую(переборщили)   
                vr = V;
            if (k == 0) // новая орбита ещё внутри
                vl = V;
            V = (vl + vr) / 2.;
        }
        W.ORB = new_A;
        W.SPEED =  V - v_p();    
    }
    return W;
}


//
double manevr(orbit A, double eps){ 
    sputnik B(A, g); 
    std::array<double,3> a = peresech(m_orb, A); // направляющий вектор пересечения орбит
    double SUM = 0.;
    if(a[0] == 0 && a[1] == 0 && a[2] == 0) // т.е. лежит в плоскости
    {
        std::cout<< "плоскость целевого эллипса совпала с настоящей" << std::endl;

        std::vector<double> x = intersection(m_orb, A, eps); // экстремальные точке
        std::vector<double> PHI; // вектор углов из которых мы видим точку пересечения
        std::vector<double> UGOL; // угол между скоростями(начальной и целевой)

        if (x.size() == 1){ // если экстремальные точки 
            std::cout<< "Добавить скорость в угле " << x[0] << "\n";
            std::cout<< "dV = " << delta_V(A, x[0])  << "km/s" << "\n";
            std::cout<< "Скорость под углом" << 0 << "(касание)" << "\n";
            SUM += abs(delta_V(A, x[0]));       
        }
        if (x.size() == 2){  // добавление двух элементов 
            int l = 0;
            if (delta_V(A, x[0]) < delta_V(A, x[1]))
                l = 0;
            else l = 1;
            std::cout<< "Добавить скорость в угле " << x[l] << "\n";
            std::cout<< "dV = " << delta_V(A, x[l])  << "km/s" << "\n";
            std::cout<< "Скорость под углом = " << acos(cos_beta(m_orb,A,x[l])) << "\n";
            SUM += abs(delta_V(A, x[l]));            
        }
        if (x.size() == 0){
            Orbit_and_V W = v_min_p(A,eps);
            std::vector<double> y = intersection(W.ORB,A,eps);
            sputnik SP(W.ORB, g);
            SUM += abs(W.SPEED);
            if(W.SPEED > 0){ // целевая орбита снаружи
                std::cout<< "Добавить скорость в перицентре dV = " << W.SPEED << "km/s" << "\n";
                std::cout<< "Получили новую орбиту с параметрами:\n" << W.ORB << "\n";
                std::cout<< "Новая орбита касается с целевой в угле: " << y[0]  << "\n";
                std::cout<< "В этой точке добавим dV = " << SP.delta_Vr(A, y[0])  << "km/s" << "\n";
                std::cout<< "Скорость под углом: " << 0 << " (касание)" << "\n";
                SUM += abs(SP.delta_V(A, y[0]));      
            }
            else // нужно уменьшать скорость
            {

                std::cout<< "Добавить скорость в перицентре dV = " << W.SPEED << "km/s" << "\n";                    
                std::cout<< "Получили новую орбиту с параметрами:\n" << W.ORB << "\n";
                
                std::cout<< "Новая орбита касается с целевой в угле: " << y[0]  << "\n";
                std::cout<< "В данной точке добавим dV = " << SP.delta_Vr(A, y[0])  << "km/s" << "\n";
                std::cout<< "Скорость под углом: " << 0 << " (касание)" << "\n";
                SUM += abs(SP.delta_Vr(A, y[0]));                      
            }
            
            std::cout<< "Перешли на целевую орбиту" << std::endl;
        }
        return SUM;        
    }

    return -1;
    }

};




//написал сужение отбиты до данной при больше окружности

#endif 