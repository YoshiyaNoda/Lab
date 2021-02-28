#include <iostream>
#include <fstream>
#include <cmath>

// 定数
constexpr double m = 0.7;
constexpr double k = 0.2;
constexpr double g = 9.81;
constexpr double dt = 0.0001;
constexpr double N = 20000; // ループ回数

// 初期条件
constexpr double t0 = 0;
constexpr double v0 = 0;

double nextApproximateV(double v) {
    return v + dt * (g - k / m * (v + dt / 2 * (g - k / m * v)));
}
double nextApproximateT(double t) {
    return t + dt;
}
double culcTheoreticalValue(double t) {
    return v0 + m * g / k *(1 - std::exp(-k / m * t));
}

int main() {
    std::ofstream outputfile("approximate.txt");
    std::ofstream tOutputfile("theoretical.txt");
    std::ofstream eOutputfile("runge_2_error.txt");
    double v = v0;
    double t = t0;
    for(int i = 0; i < N; i++) {
        outputfile << t << "\t" << v << "\n";
        tOutputfile << t << "\t" << culcTheoreticalValue(v) << "\n";
        eOutputfile << t << "\t" << (culcTheoreticalValue(v) - v) << "\n";
        v = nextApproximateV(v);
        t = nextApproximateT(t);
    }
    outputfile.close();
    tOutputfile.close();
    eOutputfile.close();
    return 0;
}
