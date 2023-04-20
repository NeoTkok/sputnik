#ifndef ORBIT_HPP
#define ORBIT_HPP

#include<cmath>

struct orbit{
    const double a; // большая полуось
    const double e; // эксцентриситет
    const double OMEGA; //долгота восходящего узла
    const double omega; // аргумент перигея
    const double i; // наклонение

    //инициализация структуры
    orbit(const double a, const double e, const double OMEGA, const double omega,
          const double i) : a(a), e(e), OMEGA(OMEGA), omega(omega), i(i) {};

    // получение радиуса перицентра
    double getRp() const{
        return a * (1 - e);
    }
    
    // получение радиуса апоцентра
    double getRa() const{
        return a * (1 + e);        
    }

    // получение фокального па
    double getP() const{
        return a * (1 - e * e);
    }

    // получение нормированных коэффициентов в уравнеии плоскоти Ax+By+Cz=0
    // D = 0, эллипс задаёт плоскость
    std::array<double, 3> getPlane() const{
        std::array<double, 3> r;
        r[0] = sin(i) * sin(OMEGA);
        r[1] = -sin(i) * cos(OMEGA);
        r[2] = cos(i);
        return r;
    }
};
const int N = 1e5;


// значение разности двух орбит в плоскости при угле фи (D - целевая орбита)
// если отрицательное, то значение целевой орбиты в данном угле ближе к планете
double delta_f(orbit D0, orbit D, double phi){
    return D0.getP() / (1 + D0.e * cos(phi)) - D.getP() / (1 + D.e * cos(phi + D0.omega - D.omega));
} 

// делит угол на .. частей и определяет части, на которых функция меняет знак 
// и отделяет эти отрезки -1, если точка является нулём, то только ей 
std::vector<double> null_f(orbit D0, orbit D){
    std::vector<double> x;
    double d = delta_f(D0,D,0);
    
    for(double alpha = 2 * M_PI / N; alpha < 2 * M_PI; alpha += 2 * M_PI / N){
        double y = delta_f(D0,D, alpha);

        if (d * y <= 0){ // условие достяжения орбиты
            if(y == 0){ // вероятность в реальной жизни такого 0, но сразу отметим 
                x.push_back(-1);
                x.push_back(alpha);
                x.push_back(-1);
            }
            // между этими значениями угла функция меняет знак
            x.push_back(-1);
            x.push_back(alpha - 2*M_PI / N);
            x.push_back(alpha);
            x.push_back(-1);
        }
        d = y;
    }
    return x;
}

// находит нули между двумя углами с точностью eps(функция в которых имеет разные знаки) 
// алгоритм основывается на методе половинного деления
double exact_null_f(const orbit D0, const orbit D, double a1, double a2, double eps){
    double l = (a1 + a2) / 2;
    double Z = delta_f(D0, D, l);
    while(Z * Z > eps * eps){
        if (Z * delta_f(D0, D, a1) > 0)
            a1 = l;
        else
            a2 = l; 
        l = (a1 + a2) / 2;
        Z = delta_f(D0, D, l);
    }
    return l;
}

/*********************************/

// аналогично но только для произвдной
double delta_df(orbit D0, orbit D, double phi){
    double y2 = 1 + D.e * cos(phi + D0.omega - D.omega);
    double y1 = 1 + D0.e * cos(phi);
    return D0.e * D0.getP() * sin(phi) / (y1*y1) - D.e * D.getP() * sin(phi + D0.omega - D.omega) / (y2 * y2);
}



std::vector<double> null_df(orbit D0, orbit D){
    std::vector<double> x;
    double d = delta_df(D0,D,0);
    
    for(double alpha = 2 * M_PI / N; alpha < 2 * M_PI; alpha += 2 * M_PI / N){
        double y = delta_df(D0,D, alpha);

        if (d * y <= 0){ // условие достяжения орбиты
            if(y == 0){ // вероятность в реальной жизни такого 0, но сразу отметим 
                x.push_back(-1);
                x.push_back(alpha);
                x.push_back(-1);
            }
            // между этими значениями угла функция меняет знак
            x.push_back(-1);
            x.push_back(alpha - 2*M_PI / N);
            x.push_back(alpha);
            x.push_back(-1);
        }
        d = y;
    }
    return x;

} 

double exact_null_df(const orbit D0, const orbit D, double a1, double a2, double eps){
    double l = (a1 + a2) / 2;
    double Z = delta_df(D0, D, l);
    while(Z * Z > eps * eps){
        if (Z * delta_df(D0, D, a1) > 0)
            a1 = l;
        else
            a2 = l; 
        l = (a1 + a2) / 2;
        Z = delta_df(D0, D, l);
    }
    return l;
}


// выдаёт нормированный вектор линии пересечения двух плоскостей(эллипс задаёт плоскость)
std::array<double, 3> peresech(const orbit A, const orbit B){
    std::array<double, 3> An = A.getPlane();
    std::array<double, 3> Bn = B.getPlane();
    std::array<double, 3> X;
    X[0] = An[1]*Bn[2] - An[2]*Bn[1];
    X[1] = - An[0]*Bn[2] + An[2]*Bn[0];
    X[2] = An[0]*Bn[1] - An[1]*Bn[0];
    double r = sqrt(X[0] * X[0] +  X[1] * X[1] + X[2] * X[2]);
    for(int i = 0; i < 3; ++i)
        X[i] /= r;
    return X;
}




#endif 
