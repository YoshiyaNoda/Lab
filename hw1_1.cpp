#include <iostream>
#include<fstream>

// 定数
constexpr double m = 0.7;
constexpr double k = 0.2;
constexpr double g = 9.81;
constexpr double dt = 0.1;

// 初期条件
constexpr double t0 = 0;
constexpr double v0 = 0;

double nextV(double v) {
    return (g - k/m * v) * dt + v;
}
double nextT(double t) {
    return t + dt;
}

int main() {
    ofstream outputfile("resule.txt");
    double v = v0;
    double t = t0;
    for(int i = 0; i < 20; i++) {
        std::cout << v << " " << t << std::endl;
        v = nextV(v);
        t = nextT(t);
    }
    return 0;
}
