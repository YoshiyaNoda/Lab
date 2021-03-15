#include <iostream>
#include <fstream>
#include <cmath>
#define GNUPLOT_PATH "gnuplot"
#include <vector>

constexpr int N = 256;

std::vector<std::vector<double>> configX() {
    std::vector<std::vector<double>> xs(N, std::vector<double>(3));
    return xs;
}
std::vector<std::vector<double>> configV() {
    std::vector<std::vector<double>> vs(N, std::vector<double>(3));
    return vs;
}
// ここで初期位置、初期速度、初期温度の設定
std::vector<std::vector<std::vector<double>>> config() {
    const std::vector<std::vector<double>> vs = configV();
    const std::vector<std::vector<double>> xs = configX();
    return {vs, xs};
}

int main() {
    const std::vector<std::vector<std::vector<double>>> ini = config();
    std::cout << ini[0].size() << std::endl;
}