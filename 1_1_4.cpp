#include <iostream>
#include <fstream>
#include <cmath>

// 定数
constexpr double m = 0.7;
constexpr double k = 0.2;
constexpr double g = 9.81;
constexpr double M = 30001;

// 初期条件
constexpr double t0 = 0;
constexpr double v0 = 0;

std::ofstream eOutputfile("youkaihou_error_114.txt");
std::ofstream inkaihouErrorOutputfile("inkaihou_error_114.txt");
std::ofstream runge2ErrorOutputfile("runge2_error_114.txt");
std::ofstream runge4ErrorOutputfile("runge4_error_114.txt");

double youkaihouNextApproximateV(double v, double dt) {
    return (g - k / m * v) * dt + v;
}
double inkaihouNextApproximateV(double v, double dt) {
    return (v + g * dt) / (1 + k / m * dt);
}
double runge2NextApproximateV(double v, double dt) {
    return v + dt * (g - k / m * (v + dt / 2 * (g - k / m * v)));
}
double runge4NextApproximateV(double v, double dt) {
    double result;
    {
        double v1 = v + dt / 2 * (g - k / m * v);
        double v2 = v + dt / 2 * (g - k / m * v1);
        double v3 = v + dt * (g - k / m * v2);
        result = v + dt / 6 * ((g - k / m * v) + 2 * (g - k / m * v1) + 2 * (g - k / m * v2) + (g - k / m * v3));
    }
    return result;
}
double nextApproximateT(double t, double dt) {
    return t + dt;
}
double culcTheoreticalValue(double t) {
    return v0 + m * g / k * (1 - std::exp(-k / m * t));
}

void youkaihou(double dt, int N) {
    double t = t0;
    double v = v0;
    for(int i = 0; i < N + 1; i++) {
        // t = 3 s
        if(i == N) {
            eOutputfile << dt << "\t" << std::abs(culcTheoreticalValue(t) - v) << "\n";
        }
        t = nextApproximateT(t, dt);
        v = youkaihouNextApproximateV(v, dt);
    }
}
void inkaihou(double dt, int N) {
    double t = t0;
    double v = v0;
    for(int i = 0; i < N + 1; i++) {
        // t = 3 s
        if(i == N) {
            inkaihouErrorOutputfile << dt << "\t" << std::abs(culcTheoreticalValue(t) - v) << "\n";
        }
        t = nextApproximateT(t, dt);
        v = inkaihouNextApproximateV(v, dt);
    }
}
void runge2(double dt, int N) {
    double t = t0;
    double v = v0;
    for(int i = 0; i < N + 1; i++) {
        // t = 3 s
        if(i == N) {
            runge2ErrorOutputfile << dt << "\t" << std::abs(culcTheoreticalValue(t) - v) << "\n";
        }
        t = nextApproximateT(t, dt);
        v = runge2NextApproximateV(v, dt);
    }
}
void runge4(double dt, int N) {
    double t = t0;
    double v = v0;
    for(int i = 0; i < N + 1; i++) {
        // t = 3 s
        if(i == N) {
            runge4ErrorOutputfile << dt << "\t" << std::abs(culcTheoreticalValue(t) - v) << "\n";
        }
        t = nextApproximateT(t, dt);
        v = runge4NextApproximateV(v, dt);
    }
}

int main() {
    for(int i = M; i > 0; i--) {
        // 切り捨てられないように
        double div = i;
        double dt = 3 / div;
        youkaihou(dt, i);
        inkaihou(dt, i);
        runge2(dt, i);
        runge4(dt, i);
    }
    eOutputfile.close();
    inkaihouErrorOutputfile.close();
    runge2ErrorOutputfile.close();
    runge4ErrorOutputfile.close();
}
