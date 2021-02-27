#include <iostream>
#include <cassert>

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

int main() {
    double v = v0;
    for(int i = 0; i < 20; i++) {
        v = nextV(v);
        std::cout << v << std::endl;
    }
    return 0;
}
