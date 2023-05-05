#ifndef SPUTNIK_HPP
#define SPUTNIK_HPP

#include "Orbit.hpp"

class sputnik {
private:
    orbit m_orb;
    double g;
public:
    sputnik(const orbit& A, const double G) : m_orb(A), g(G) {};

// получение энергии системы
    double get_h() const{
        return - g * (1-m_orb.e * m_orb.e) / m_orb.p;
    }

// радиальная скорость в зависимости от угла(от перицентра против ч.с.)
    double v_r(const double alpha) const{
        return sqrt(g/m_orb.p) * m_orb.e * sin(alpha);
    }

// трансверсальная скорость ....-//-
    double v_tau(const double alpha) const{
        return sqrt(g / m_orb.p) * (1 + m_orb.e * cos(alpha));
    }

//***********
// скорость в перицентре(просто удобно отдельно так сделать)
    double v_p() const{
        return (1 + m_orb.e) * sqrt(g / m_orb.p);
    }

//***************
// новый фокальный параметр при изменении скорости в перицентре на dv
    double new_p_p(const double dv) const{
        double c = (v_p() + dv) * m_orb.p / (1. + m_orb.e);
        return c * c / g;
    }

//****************
// новый e при изменении скорости в перицентре на dv
    double new_e_p(double dv) const{
        return (v_p() + dv) * (v_p() + dv) * m_orb.p / (1. + m_orb.e) / g - 1.;
    }

// новый фокальный параметр при изменении скорости в апоцентре на dv
    double new_p_a(const double dv) const{
        double c = (v_p() + dv) * m_orb.p / (1. - m_orb.e);
        return c * c / g;
    }
//****************
// новый e при изменении скорости в апоцентре на dv
    double new_e_a(const double dv) const{
        return 1 - (v_p() + dv) * (v_p() + dv) * m_orb.p / (1. - m_orb.e) / g ;
    }


// получение скорости при заданном угле
    double get_v(const double phi) const{
        return sqrt(get_h() + 2. * g / m_orb.getR(phi)) ;
    }


// получение характерной скорости в данном угле начальной орбиты 
// при заданном угле пересечения с целевой
    double delta_V(const orbit& A, const double phi) const{
        sputnik B(A,g);
        double b = cos_beta(m_orb, A, phi); // косинус угола межнду начальной скоростью и конечой
        double v0 = get_v(phi); // как раз таки начальная скорость
        double v = B.get_v(phi); // скорость, которую будет иметь КА на целевой орбите

        return sqrt(v*v + v0*v0 - 2*v*v0*b); // т. Косинусов
    }

// получение скоростей при ОЧЕНЬ малых углах(где касательные)
    double delta_Vr(const orbit& A, const double phi) const{
        sputnik B(A,g);
        double v0 = get_v(phi); 
        double v = B.get_v(phi); 
        double sign = (m_orb.S * A.S > 0.)?1.:-1.;
        return v - sign*v0;
    }





// а вот эта функция...
// происходит ситуация, когда на плоскости орбиты не пересекаются
// (например целевая орбита находится внутри или снаружи)
// поэтом данная функция высчитывает скороть которую надо добавить или убавить в
// перицентре чтобы произошло касание
// (здесь есть несколько подводных камней, которые стоит обсудить лично 
// не только что это не прям оптимально, но и при переходе на очень отдаленную орбиту
// (которая около в 6-7 раз дальше) перейти не получится, т к там уже выходим на параболу и тд...)
Orbit_and_V v_min_p(const orbit& A, const double eps) const{
    Orbit_and_V W(m_orb,0.); // эту сущность я собираюсь возврящать
    if(delta_f(m_orb,A,0.) > 0.){ // целевая орбита внутри
        // далее мне надо проверить, что данная орбита пересекает окружность
        // радиуса перецентра или нет, так как в противном случае перицентр и апоцентр
        // меняются местами и при изменении скорости е и р меняются уже по-другому
        // для этого создаю эту окружность и смотрю колличесво пересечений:
        double rp = m_orb.getRp(); 
        double Vo = sqrt(g / rp);
        orbit new_A = m_orb;
        new_A.e = 0;
        new_A.p = rp;
        int n = intersection(new_A, A, eps).size(); // колличество пересечений
        if (n == 1) // +целевая орбита касается окружности радиуса перицентр    
        {
            W.ORB = new_A; // возвращаем орбиту - окружность
            W.SPEED =  Vo - v_p(); // а скорость - изменение затраченное на маневр
        }
        if (n == 0){ // целевая орбита внутри окружности разности перецентр 
            double vl = 0.; // левая граница скорости
            double vr = Vo; // правая граница скорости 
            double V = (vl + vr) / 2.; // серединка скорости
            
            new_A.omega += M_PI; // вотздесь требуется повернуть главную ось.. 
            if (new_A.omega >= 2*M_PI)
                new_A.omega -= 2 * M_PI;
            sputnik G(new_A, g);
            orbit new_B = new_A; // сейчас находимся на круговой орбите            
            while(true){ 
                // стоит заметить, что скорость добавляем в апоцентре                    
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
            return W; // всё)
        }
        if (n == 2){  // пересекает окружность    
            double vl = Vo; // левая граница скорости
            double vr = v_p(); // правая граница скорости
            double V = (vl + vr) / 2.; // серединка скорости 
            
            orbit new_B = m_orb; // новая орбита после добавления скорости
            while(true){
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
        // правую границу сделаем близ лежащую к параболе (достигается при числе не 1.87 а 2.)
        double vr = sqrt((1.87)*g*(1+m_orb.e)/m_orb.p); // при более больших чем 1.87-1.9 ломается
        // т к отбита становится сильно вытянутой и требуется уже другие маневры
            
        new_A.e = new_e_p(vr - v_p()); // новый эксцентриситет
        new_A.p = new_p_p(vr - v_p()); // новый фокальный параметр            
        int k = intersection(new_A, A, eps).size();
            
        if (k == 0){
            std::cout << "нельзя добавить в перцентр скорость" << std::endl;
            return W; // так не идеально сделать, но зато легко отловить
        }
        double V = (vl + vr) / 2.; // серединка скорости 

        while(true){ // 2^15 ошибка убывает
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

// функция вернёт ноую орбиту и новую скорость
Orbit_and_V povorot(const orbit& A, const double eps) const{
    Orbit_and_V W(m_orb,0);
    std::array<double, 3> a = peresech(m_orb,A, eps);
    double phi_1 = PHI(a,m_orb,eps); // угол на вектор пересечения у начальной орбиты
    double phi_2 = PHI(a,A,eps);  // -//- у целевой орбиты
   
    std::array<double,3>  x = m_orb.a_r(phi_1);

    W.ORB.i = A.i; // после поворота произойдет так
    // а теперь сам поворот
    sputnik S(A,g);
    double Vo = get_v(phi_1);
    double V = S.get_v(phi_2);
    double COS = cos_gama(m_orb,A,eps);
    W.SPEED = sqrt(V*V + Vo*Vo - 2*Vo*V*COS);
    return W;
}






// это уже основная функция маневра
double manevr(const orbit& A, const double eps) const{

    std::array<double,3> a = peresech(m_orb, A, eps); // направляющий вектор пересечения орбит
    double SUM = 0.; // суммарня характерная скорость(которую требуется найти)     
    Orbit_and_V H = povorot(A,eps);
    orbit Z = m_orb; 

    if(a[0] != 0 or a[1] != 0 or a[2] != 0) // 
    {
        Z = H.ORB;
        std::cout << "Плоскости орбит не совпадают!" << std::endl;
        std::cout << "Поддадим скорость в угле: " << PHI(a,m_orb,eps)<< " начальной орбиты" << std::endl; 
        std::cout << "Под углом: " << acos(cos_gama(m_orb,A,eps)) <<  " к целевой орбите" << std::endl;
        SUM += H.SPEED;
        std::cout << "Характерная скорость данного поворота: " << H.SPEED << "m/s" << std::endl;
        std::cout << "Таким образом мы перешли на другую орбиту, лежащу в плоскости целевой орбиты:" << std::endl << H.ORB;
    }

    sputnik B(Z, g);
    { // маневр происходит теперь вплоскости
        std::cout<< "Плоскость целевого эллипса совпала с настоящей!" << std::endl;
        std::vector<double> x = intersection(Z, A, eps); // экстремальные точке
        // вот здесь 3 варианта: касание, бесконечно много(совпадение); 0 и 2 пересечений 
        if (x.size() == 1){ // касание 
            std::cout<< "Добавить скорость в угле " << x[0] << "\n";
            std::cout<< "dV = " << B.delta_Vr(A, x[0])  << "m/s" << "\n";
            std::cout<< "Скорость под углом: " << 0 << "(касание)" << "\n";
            SUM += abs(B.delta_V(A, x[0]));       
        }
        if (x.size() == 2){ // орбиты пересекаются в 2х точках 
            int l = 0;
            if (B.delta_V(A, x[0]) < B.delta_V(A, x[1])) // из двух точек смотрим более оптимальную
                l = 0;
            else l = 1;
            std::cout<< "Добавить скорость в угле: " << x[l] << "\n";
            std::cout<< "dV = " << B.delta_V(A, x[l])  << "m/s" << "\n";
            std::cout<< "Скорость под углом = " << acos(cos_beta(Z,A,x[l])) << "\n";
            SUM += abs(B.delta_V(A, x[l]));            
        }
        if (x.size() == 0){ // орбиты не пересекаются
            Orbit_and_V W = B.v_min_p(A,eps);
            std::vector<double> y = intersection(W.ORB,A,eps);
            sputnik SP(W.ORB, g);
            SUM += abs(W.SPEED);
            if(W.SPEED > 0){ // целевая орбита снаружи
                std::cout<< "Добавить скорость в перицентре dV = " << W.SPEED << "km/s" << "\n";
                std::cout<< "Получили новую орбиту с параметрами:\n" << W.ORB << "\n";
                std::cout<< "Новая орбита касается с целевой в угле: " << y[0]  << "\n";
                std::cout<< "В этой точке добавим dV = " << SP.delta_Vr(A, y[0])  << "km/s" << "\n";
                std::cout<< "Скорость под углом: " << 0 << " (касание)" << "\n";
                SUM += abs(SP.delta_Vr(A, y[0]));      
            }
            if(W.SPEED < 0) // нужно уменьшать скорость (целевая внутри)
            {

                std::cout<< "Добавить скорость в перицентре dV = " << W.SPEED << "m/s" << "\n";                    
                std::cout<< "Получили новую орбиту с параметрами:\n" << W.ORB << "\n";
                
                std::cout<< "Новая орбита касается с целевой в угле: " << y[0]  << "\n";
                std::cout<< "В данной точке добавим dV = " << SP.delta_Vr(A, y[0])  << "m/s" << "\n";
                std::cout<< "Скорость под углом: " << 0 << " (касание)" << "\n";
                SUM += abs(SP.delta_Vr(A, y[0]));                      
            }
            
            std::cout<< "Перешли на целевую орбиту" << std::endl;
        }
        return SUM;        
    }

    return SUM;
    }

};




//написал сужение отбиты до данной при больше окружности

#endif 