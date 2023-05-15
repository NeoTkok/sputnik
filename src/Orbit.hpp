#ifndef ORBIT_HPP
#define ORBIT_HPP

#include<cmath>

struct orbit{
    double p; // параметр орбиты
    double e; // эксцентриситет
    double OMEGA; //долгота восходящего узла
    double omega; // аргумент перигея
    double i; // наклонение
    std::array<double, 3>  S = {1., 1. ,1.}; // "положиельная ориентация" движени КА вдоль орбиты
    //инициализация структуры
    orbit(const double p, const double e, const double OMEGA, const double omega,
          const double i) : p(p), e(e), OMEGA(OMEGA), omega(omega), i(i) {};

    orbit(const double p, const double e, const double OMEGA, const double omega, const double i,
          const std::array<double,3> S) : p(p), e(e), OMEGA(OMEGA), omega(omega), i(i), S(S) {};


    // получение радиуса перицентра
    double getRp() const{
        return p / (1 + e);
    }
    
    // получение радиуса апоцентра
    double getRa() const{
        return p / (1 - e);        
    }

    double getR(const double phi) const{
        return p / (1 + e * cos(phi));
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
    // направление вектора радиального направления в пространстве   
    std::array<double, 3> a_r(const double phi) const{
        std::array<double, 3> a = {0,0,0};
        a[0] = getR(phi)*(cos(OMEGA)*cos(omega+phi)-sin(OMEGA)*sin(omega+phi)*cos(i));
        a[1] = getR(phi)*(sin(OMEGA)*cos(omega+phi)+cos(OMEGA)*sin(omega+phi)*cos(i));
        a[2] = getR(phi)*sin(omega+phi)*sin(i);    
        return a;
    }
};
const double N = 1000;


// разность двух векторов в пространстве 
double operator*(const std::array<double,3>& A, const std::array<double,3>& B){
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}
// выдаёт нормированный вектор линии пересечения двух плоскостей(эллипс задаёт плоскость)
std::array<double, 3> peresech(const orbit& A, const orbit& B, const double eps){
    std::array<double, 3> An = A.getPlane(); // коэффициенты 1 плоскости
    std::array<double, 3> Bn = B.getPlane(); // коэффициенты 2 плоскости
    std::array<double, 3> X; // сюда я буду складывать направляющие векторы пересечения
    X[0] = An[1]*Bn[2] - An[2]*Bn[1];   // 
    X[1] = - An[0]*Bn[2] + An[2]*Bn[0]; //здесь появеляется векторное произведение 
    X[2] = An[0]*Bn[1] - An[1]*Bn[0];   //
    double r = sqrt(X[0] * X[0] +  X[1] * X[1] + X[2] * X[2]);
    if(r < eps)
        {
            X[0] = 0.; X[1] = 0.; X[2] = 0.;
            return X;
        }       
    for(int i = 0; i < 3; ++i)
        X[i] /= r;
    return X;
}



// значение разности двух орбит в плоскости при угле фи (D - целевая орбита)
// если отрицательное, то значение целевой орбиты в данном угле ближе к планете
double delta_f(const orbit& D0, const orbit& D, const double phi){
    return D0.getR(phi) - D.getR(phi + D0.omega - D.omega);
} 

// аналогично но только для произвдной
double delta_df(const orbit& D0, const orbit& D, const double phi){
    double y2 = 1 + D.e * cos(phi + D0.omega - D.omega);
    double y1 = 1 + D0.e * cos(phi);
    return D0.e * D0.p * sin(phi) / (y1 * y1) - D.e * D.p * sin(phi + D0.omega - D.omega) / (y2 * y2);
}


// функция возвращающая отрезки угла, где производная функции меняет знак(т.е.  = 0)
// в таком формате : (-1,х1,у1,-1,-1,х2,у2,-1)
// таких отрезков не больше чем 4, но в основном только 2))
std::vector<double> null_df(const orbit& D0, const orbit& D){
    std::vector<double> x;
    double d = delta_df(D0,D,0.); // получаем значение производной функции в перицентре
    
    for(double alpha = 2 * M_PI / N; alpha < 2 * M_PI; alpha += 2 * M_PI / N){
        double y = delta_df(D0,D, alpha);
        if (d * y <= 0){ 
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
// (метод половинного деления)
double exact_null_df(const orbit& D0, const orbit& D, double a1, double a2, const double eps){
    double l = (a1 + a2) / 2;
    double Z = delta_df(D0, D, l);
    while(abs(Z) > eps){
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
std::vector<double> extr(const orbit& D0, const orbit& D, const double eps){
    std::vector<double> n_df = null_df(D0, D); // находим отрезки где производная меняет знак
    std::vector<double> ex_n_df;
    
    for(int i = 1; i < n_df.size() - 2; ++i) // получили вектоpa значений угла в которых происходит 
        if(n_df[i - 1] + n_df[i + 2] == -2 && n_df[i] != -1 && n_df[i + 1] != 0) // -//- но только с поиском
            ex_n_df.push_back(exact_null_df(D0, D, n_df[i], n_df[i + 1], eps));
    // зачем я делаю эту вещь - тяжело описать словами
    // поэтому расскажу лично
    if (ex_n_df.size() > 2){
        std::vector<double> x(2);
        if(ex_n_df.size() > 4)
            return ex_n_df;

        x[0] = ex_n_df[0];
        x[1] = ex_n_df[2];
        return x;
    }   
    return ex_n_df;

}

// находит нули между двумя углами с точностью eps(функция в которых имеет разные знаки) 
// алгоритм основывается на методе половинного деления
double exact_null_f(const orbit& D0, const orbit& D, double a1, double a2, const double eps){
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


// функция, проверяющая является ли касанием двух орбит в некотором phi
// при условии, что phi - экстремум функции f
bool touch(const orbit& D0, const orbit& D, const double phi, const double eps){
    if (abs(delta_f(D0,D,phi)) / D0.getR(phi) < eps) 
        return 1;
    return 0;
}

// выдаёт массив размера колличества пересечений (размер 0, 1 или 2)
// если пересечений больше 2 (т е идет наложение) функция выдаёт 3 нуля
std::vector<double> intersection(const orbit& D0, const orbit& D, const double eps){
    std::vector<double> ex_n_df = extr(D0, D, eps); 
    std::vector<double> ex_n_f; // точные точки пересечения или касания

    if(ex_n_df.size() > 2){
        std::vector<double> Z = {0. ,0. ,0. }; // происходит наложение
        return Z;
    }

    if (delta_f(D0,D,ex_n_df[0]) * delta_f(D0,D,ex_n_df[1]) < 0) // меняет знак => 2 пересечения - 2 общие точки
    {
        ex_n_f.push_back(exact_null_f(D0, D, ex_n_df[0], ex_n_df[1], eps));        
        ex_n_f.push_back(exact_null_f(D0, D, ex_n_df[1], ex_n_df[0] + 2 * M_PI, eps));
        if (ex_n_f[1] >= 2*M_PI) 
            ex_n_f[1] -= 2*M_PI;
        return ex_n_f;
    }

    if (touch(D0,D,ex_n_df[0],eps) == 1){ // касание в 1 точке - только одна общая точка
        ex_n_f.push_back(ex_n_df[0]);
        return ex_n_f;
    }

    if (touch(D0,D,ex_n_df[1],eps) == 1){ // касание во второй точке - только одна общая точка 
        ex_n_f.push_back(ex_n_df[1]);
        return ex_n_f;
    }
    
    return ex_n_f;
}


// функция вычисляет косинус угла между орбитами,
// которые пересекаются под углом phi(от планеты)
// данный алгоритм дает плохую точность при очень близких углах(касание)-
// для решения этой проблемы косинус будет просто равен +-1
// а для больших углов (близ 90 градусов) хорошоая точноть,но такую ситуации в реальном мире очень сложно повторить
// тк под прямым углом никакие орбиты не пересекаются(у одной е->0 а у другой у->1 и еще пару условий)
double cos_beta(const orbit& D0, const orbit& D, const double phi){

    double Rphi = D0.getR(phi); // расстояние до точки пересечения
    double alpha = 0.0001; // угол отступа от пересечения
    double R0 = D0.getR(phi + alpha);  // расстояние до данной отбиты отклоненной от пересеч на alpha
    double R = D.getR(phi + D0.omega - D.omega + alpha); // аналогично, но только для целевой отбиты
    alpha = (1. - alpha * alpha / 2 + alpha*alpha*alpha*alpha / 24); // разложение косинуса в нуле

    double d0 = (R0 * R0 + Rphi * Rphi - 2 * R0 * Rphi * alpha); // растояние на начальной орбите 
    double d =(R * R + Rphi * Rphi - 2 * R * Rphi * alpha); // аналогично но для целевой
    
    double sign = (D0.S*D.S > 0.)?1.:-1.;

    return (R * R0 + Rphi * Rphi - sign * Rphi * (R0 + R) * alpha ) / sqrt(d * d0); //возвращаю косинус)

    // для улучшения малых улов есть одна идейка
    
}

// разность двух векторов в пространстве 
std::array<double, 3> operator-(const std::array<double,3>& A, const std::array<double,3>& B){
    std::array<double, 3> x = A;
    x[0] -= B[0];
    x[1] -= B[1];
    x[2] -= B[2];
    return x;
}

// квадрат вектора
double mod(const std::array<double,3>& A){
    return A[0]*A[0]+A[1]*A[1]+A[2]*A[2];    
}


// функция получающая максимально близкий вектор орбиты, соноправленный с а
double PHI(const std::array<double,3>& a, const orbit& A, const double eps){
    double l = 0.; // нижняя граница    
    double val_l = mod(a - A.a_r(0.)); // расстояние снизу
    double r = 2. * M_PI; // верхняя граница
    double val_r = mod(a - A.a_r(2.*M_PI)); // расстояние сверху
    //std::cout << val_r << " " << val_l << std::endl; 
    if (val_r > mod(a-A.a_r(r - eps))){  // если мы можем уменьшать верхнюю границу
    // то в силу монатонности найдем данный угол
        while(val_r > mod(a-A.a_r(r -  eps)) and r >= 0){
            r -= eps;
            val_r = mod(a-A.a_r(r));
//            std::cout << r << std::endl;

        }
        return r - A.omega;
    }
    else{ // если верхняя граница не уменьшается, то будем нижнюю увеличивать, пока можем
        while(val_l > mod(a-A.a_r(l + eps)) and l <= 2 * M_PI){
            l += eps;
            val_l = mod(a-A.a_r(l));
  //          std::cout << l << std::endl;

        }
        return l - A.omega;
    }
}

// угол между направлениями в пространстве 
double cos_gama(const orbit& D0, const orbit& D, const double eps){
    
    std::array<double, 3> a = peresech(D0,D, eps);
    double phi_1 = PHI(a,D0,eps); // угол на вектор пересечения у начальной орбиты
    double phi_2 = PHI(a,D,eps);  // -//- у целевой орбиты

    std::array<double, 3> Rphi_1 = D0.a_r(phi_1); //расстояние до первой орбиты
    std::array<double, 3> Rphi_2 = D.a_r(phi_2); // -//- до второй

//    std::cout << Rphi_1[0] << " " << Rphi_1[1] << " " << Rphi_1[2]  << std::endl; 
//    std::cout << Rphi_2[0] << " " << Rphi_2[1] << " " << Rphi_2[2]  << std::endl; 

    double alpha = 0.0001;
    std::array<double, 3> R1 = D0.a_r(phi_1 + alpha) - Rphi_1;
    std::array<double, 3> R2 = D.a_r(phi_2 + alpha) - Rphi_2;
    
    //std::cout << R1[0] << " " << R1[1] << " " << R1[2]  << std::endl; 
    //std::cout << R2[0] << " " << R2[1] << " " << R2[2]  << std::endl; 
    double sign = (D0.S*D.S > 0.)?1.:-1.;


    return sign*abs(R1[0]*R2[0] + R1[1]*R2[1] + R1[2]*R2[2])/sqrt(mod(R1) * mod(R2));
}





// наглядный вывод всех харрактеристик орбиты
std::ostream& operator<<(std::ostream& os, const orbit& o) {
    os << "############################" << std::endl;
    os << "| Параметр орбиты          | " << o.p << std::endl;
    os << "| Эксцентриситет           | " << o.e << std::endl;
    os << "| Долгота восходящего узла | " << o.OMEGA << std::endl;
    os << "| Аргумент перигея         | " << o.omega << std::endl;
    os << "| Наклонение               | " << o.i << std::endl;
    os << "############################" << std::endl;
    return os;
}



#endif 
