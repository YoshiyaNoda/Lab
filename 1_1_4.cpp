#include <iostream>
#include <fstream>
#include <cmath>

// 定数
constexpr double m = 0.7;
constexpr double k = 0.2;
constexpr double g = 9.81;

// 初期条件
constexpr double t0 = 0;
constexpr double v0 = 0;

double youkaihouNextApproximateV(double v, double dt) {
    return (g - k/m * v) * dt + v;
}
double nextApproximateT(double t, double dt) {
    return t + dt;
}
double culcTheoreticalValue(double t) {
    return v0 + m * g / k *(1 - std::exp(-k / m * t));
}

void youkaihou(double dt) {
    double t = t0;
    double v = v0;
    double t = t0;
    std::ofstream eOutputfile("youkaihou_error_114.txt");
    // t = 3 s となるようにループの回数を決める
    for(int i = 0; t < 30001; i++) {
        // t = 3 s
        if(i == 30000) {
            eOutputfile << t << "\t" << std::abs(culcTheoreticalValue(t) - v) << "\n";
        }
        v = youkaihouNextApproximateV(v, dt);
        t = nextApproximateT(t, dt);
    }
    eOutputfile.close();
}

int main() {
    double dt = 0.0001;
    for(i = 0; i < 100; i++) {
        youkaihou(dt)
        dt += 0.01;
    }
}
