#include <iostream>
#include <fstream>
#include <cmath>

// 定数
constexpr double m = 0.7;
constexpr double k = 0.2;
constexpr double g = 9.81;
constexpr double dt = 0.0001;
constexpr double N = 30001; // ループ回数

// 初期条件
constexpr double t0 = 0;
constexpr double v0 = 0;

double nextApproximateV(double v) {
    double result;
    {
        double v1 = v + dt / 2 * (g - k / m * v);
        double v2 = v + dt / 2 * (g - k / m * v1);
        double v3 = v + dt * (g - k / m * v2);
        result = v + dt / 6 * ((g - k / m * v) + 2 * (g - k / m * v1) + 2 * (g - k / m * v2) + (g - k / m * v3));
    }
    return result;
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
    std::ofstream eOutputfile("runge_4_error.txt");
    double v = v0;
    double t = t0;
    for(int i = 0; i < N; i++) {
        outputfile << t << "\t" << v << "\n";
        tOutputfile << t << "\t" << culcTheoreticalValue(t) << "\n";
        eOutputfile << t << "\t" << (culcTheoreticalValue(t) - v) << "\n";
        if(i == 30000) {
            std::cout << "t = 3.0 s\n" << culcTheoreticalValue(t) - v << std::endl;
        }
        v = nextApproximateV(v);
        t = nextApproximateT(t);
    }
    outputfile.close();
    tOutputfile.close();
    eOutputfile.close();
    return 0;
}
