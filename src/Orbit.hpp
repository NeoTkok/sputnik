#ifndef ORBIT_HPP
#define ORBIT_HPP

#include<cmath>

struct orbit{
    double a; // большая полуось
    double e; // эксцентриситет
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

    double getR(const double phi){
        return getP() / (1 + e * cos(phi));
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
const int N = 20;



/*********************************/

// аналогично но только для произвдной
double delta_df(orbit D0, orbit D, double phi){
    double y2 = 1 + D.e * cos(phi + D0.omega - D.omega);
    double y1 = 1 + D0.e * cos(phi);
    return D0.e * D0.getP() * sin(phi) / (y1*y1) - D.e * D.getP() * sin(phi + D0.omega - D.omega) / (y2 * y2);
}


//функция возвращающая отрезки угла, где производная функции  = 0
// в таком формате : (-1,х1,у1,-1,-1,х2,у2,-1)
std::vector<double> null_df(orbit D0, orbit D){
    std::vector<double> x;
    double d = delta_df(D0,D,0);
    
    for(double alpha = 2 * M_PI / N; alpha < 2 * M_PI; alpha += 2 * M_PI / N){
        double y = delta_df(D0,D, alpha);
        if (d * y <= 0){ // условие достяжения орбиты
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

// в данном отрезке находит точный нуль производной
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

//фукция, выдающая вектор экстремальных точек с ошибкой eps 
std::vector<double> extr(orbit D0, orbit D, double eps){
    std::vector<double> n_df = null_df(D0, D); // находим отрезки где производная меняет знак
    std::vector<double> ex_n_df;
    for(int i = 1; i < n_df.size() - 2; ++i) // получили вектоpa значений угла в которых происходит 
        if(n_df[i - 1] + n_df[i + 2] == -2 && n_df[i] != -1 && n_df[i + 1] != 0) // -//- но только с поиском
            ex_n_df.push_back(exact_null_df(D0, D, n_df[i], n_df[i + 1], eps));
    return ex_n_df;
}





// значение разности двух орбит в плоскости при угле фи (D - целевая орбита)
// если отрицательное, то значение целевой орбиты в данном угле ближе к планете
double delta_f(orbit D0, orbit D, double phi){
    return D0.getR(phi) - D.getR(phi + D0.omega - D.omega);
} 

// находит нули между двумя углами с точностью eps(функция в которых имеет разные знаки) 
// алгоритм основывается на методе половинного деления
double exact_null_f(const orbit D0, const orbit D, double a1, double a2, double eps){
    double l = (a1 + a2) / 2;
    double Z = delta_f(D0, D, l);
    while(abs(Z) > eps){
        if (Z * delta_f(D0, D, a1) > 0)
            a1 = l;
        else
            a2 = l; 
        l = (a1 + a2) / 2;
        Z = delta_f(D0, D, l);
    }
    return l;
}


// функция, проверяющая является ли касанием двух орбит
// при условии, что phi - экстремум функции f
bool touch(orbit D0, orbit D, double phi, double eps){
    if (abs(delta_f(D0,D,phi)) / D0.getR(phi) < eps) 
        return 1;
    return 0;
}

// выдаёт массив размера колличества пересечений (либо размер 0, 1 , 2, ...)
// если пересечений больше 2(т е идет наложение) функция выдаёт 3 нуля
std::vector<double> intersection(orbit D0, orbit D, double eps){
    std::vector<double> ex_n_df = extr(D0, D, eps); // (размер 2 если не накладывается)
    if(ex_n_df.size() > 2){
        std::vector<double> Z = {0. ,0. ,0. }; // происходит наложение
        return Z;
    }

    std::vector<double> ex_n_f; // точные точки пересечения или касания
    if (touch(D0,D,ex_n_df[0],eps*100) == 1){ // касание в 1 точке - только одна общая точка
        ex_n_f.push_back(ex_n_df[0]);
        return ex_n_f;
    }

    if (touch(D0,D,ex_n_df[1],eps*100) == 1){ // касание во второй точке - только одна общая точка 
        ex_n_f.push_back(ex_n_df[1]);
        return ex_n_f;
    }

    if (delta_f(D0,D,ex_n_df[0]) * delta_f(D0,D,ex_n_df[1]) < 0) // меняет знак => 2 пересечения - 2 общие точки
    {
        ex_n_f.push_back(exact_null_f(D0, D, ex_n_df[0], ex_n_df[1], eps));        
        ex_n_f.push_back(exact_null_f(D0, D, ex_n_df[1], ex_n_df[0] + 2 * M_PI, eps));
        if (ex_n_f[1] >= 2*M_PI) 
            ex_n_f[1] -= 2*M_PI;
        return ex_n_f;
    }
    
    return ex_n_f;
}

// функция вычисляет косинус угла между орбитами, которые пересекаются под углом phi(от планеты)
double cos_beta(orbit D0, orbit D, double phi){

    double Rphi = D0.getR(phi); // расстояние до точки пересечения
    double alpha = 0.0001; // угол отступа от пересечения
    double R0 = D0.getR(phi + alpha);  // расстояние до данной отбиты отклоненной от пересеч на alpha
    double R = D.getR(phi + D0.omega - D.omega + alpha); // аналогично, но только для целевой отбиты
    alpha = (1. - alpha * alpha / 2 + alpha*alpha*alpha*alpha / 24); // разложение косинуса в нуле

    double d0 = (R0 * R0 + Rphi * Rphi - 2 * R0 * Rphi * alpha); // растояние на начальной орбите 
    double d =(R * R + Rphi * Rphi - 2 * R * Rphi * alpha); // аналогично но для целевой
    
    return (R * R0 + Rphi * Rphi - Rphi * (R0 + R) * alpha ) / sqrt(d * d0); //возвращаю косинус)

    // для улучшения малых улов есть одна идейка
    
}

// выдаёт нормированный вектор линии пересечения двух плоскостей(эллипс задаёт плоскость)
std::array<double, 3> peresech(const orbit A, const orbit B){
    std::array<double, 3> An = A.getPlane(); // коэффициенты 1 плоскости
    std::array<double, 3> Bn = B.getPlane(); // коэффициенты 2 плоскости
    std::array<double, 3> X; // сюда я буду складывать направляющие векторы пересечения
    X[0] = An[1]*Bn[2] - An[2]*Bn[1];   // 
    X[1] = - An[0]*Bn[2] + An[2]*Bn[0]; //здесь появеляется векторное произведение 
    X[2] = An[0]*Bn[1] - An[1]*Bn[0];   //
    double r = sqrt(X[0] * X[0] +  X[1] * X[1] + X[2] * X[2]);
    if(r <= 1e-15)
        {
            X[0] = 0.; X[1] = 0.; X[2] = 0.;
            return X;
        }       
    for(int i = 0; i < 3; ++i)
        X[i] /= r;
    return X;
}


#endif 
