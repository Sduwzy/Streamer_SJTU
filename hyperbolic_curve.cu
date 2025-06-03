#include "hyperbolic_curve.h"
#include <cmath>    // 若需 sqrt, pow, log, etc. 可以保留

double hyperbolic_curve(double Vele, double cr, double b, double r, double z)
{
    // (z/b)^2 - (r/a)^2 = 1
    // Q : 焦点
    // (r, z)を通る, 焦点Qの双曲線のパラメータ
    // thetaI : 电极, theta : 求めたい点
    double a, a2, b2, Q;
    double thetaI, theta;
    double V;

    a = std::sqrt(cr * b);

    // 电极内( (z/b)^2 - (r/a)^2 >= 1 ) 则返回电极电压
    if (std::pow(z/b, 2) - std::pow(r/a, 2) >= 1) {
        V = Vele;
    }
    // 若 z == 0 => 接地(0)
    else if (z == 0) {
        V = 0.0;
    }
    else {
        // 计算
        thetaI = std::atan2(a, b); // 电极处的角度
        Q = std::sqrt(b*b + a*a);  // 焦点

        // 计算 a2
        // a2 = [ -(r^2 + z^2 - Q^2) + sqrt( (r^2 + z^2 - Q^2)^2 - 4*(-r^2 * Q^2) ) ] / 2
        double inside = (-(r*r + z*z - Q*Q) 
                         + std::sqrt(std::pow((r*r + z*z - Q*Q), 2) 
                           - 4*(-r*r * Q*Q))) / 2.0;
        a2 = std::sqrt(inside);

        // b2 = sqrt(Q^2 - a2^2)
        b2 = std::sqrt(Q*Q - a2*a2);

        // theta = atan2(a2, b2)
        theta = std::atan2(a2, b2);

        // V = Vele * log( 1 / tan(theta/2) ) / log( 1 / tan(thetaI/2) )
        double numerator   = std::log(1.0 / std::tan(theta  / 2.0));
        double denominator = std::log(1.0 / std::tan(thetaI / 2.0));
        V = Vele * (numerator / denominator);
    }

    return V;
}
