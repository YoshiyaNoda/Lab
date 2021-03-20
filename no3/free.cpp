#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#define GNUPLOT_PATH "gnuplot"
#include <vector>

constexpr int N = 256;
constexpr double a = 0.00001;
constexpr int rand_min = 100;
constexpr int rand_max = -100;

std::random_device rd;
std::default_random_engine eng(rd());
std::uniform_real_distribution<double> distr(rand_min, rand_max);

std::vector<double> vRamdomConfig() {
    return {distr(eng), distr(eng), distr(eng)};
}
std::vector<std::vector<std::vector<double>>> culcHexData(std::vector<double> head_loc) {
    const std::vector<std::vector<double>> vs = {vRamdomConfig(), vRamdomConfig(), vRamdomConfig(), vRamdomConfig()};
    const std::vector<double> r1 = {head_loc[0] + a, head_loc[1] + a, head_loc[2]};
    const std::vector<double> r2 = {head_loc[0], head_loc[1] + a, head_loc[2] + a};
    const std::vector<double> r3 = {head_loc[0] + a, head_loc[1], head_loc[2] + a};
    const std::vector<std::vector<std::vector<double>>> data = {{head_loc, vs[0]}, {r1, vs[1]}, {r2, vs[2]}, {r3, vs[3]}};
    return data;
}
// ここで初期位置、初期速度、初期温度の設定
std::vector<std::vector<std::vector<double>>> config() {
    std::vector<std::vector<std::vector<double>>> molecular_data(N, std::vector<std::vector<double>>(2));
    // どの分子か->座標か速度か->その値(3次元)
    std::vector<double> head_loc = {0, 0, 0};
    // ここでは正四面体の頭の座標を更新していいく 2*aずつヅラせばいい
    for(int i = 0; i < N / 4; i++) {
        const std::vector<std::vector<std::vector<double>>> hex_data = culcHexData(head_loc);
        molecular_data[i * 4] = hex_data[0];
        molecular_data[i * 4 + 1] = hex_data[1];
        molecular_data[i * 4 + 2] = hex_data[2];
        molecular_data[i * 4 + 3] = hex_data[3];
        if(i % 4 == 3) {
            head_loc[0] = 0;
            head_loc[1] += 2 * a;
        } else {
            head_loc[0] += 2 * a;
        }
        if(i % 16 == 15) {
            head_loc[1] = 0;
            head_loc[2] += 2 * a;
        }
    }

    return molecular_data;
}

int main() {
    const std::vector<std::vector<std::vector<double>>> ini = config();
    std::cout << ini.size() << std::endl;
}