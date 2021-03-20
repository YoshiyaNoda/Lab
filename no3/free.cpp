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
std::vector<std::vector<double>> culcHexLoc(std::vector<double> head_loc) {
    const std::vector<double> r1 = {head_loc[0] + a, head_loc[1] + a, head_loc[2]};
    const std::vector<double> r2 = {head_loc[0], head_loc[1] + a, head_loc[2] + a};
    const std::vector<double> r3 = {head_loc[0] + a, head_loc[1], head_loc[2] + a};
    return {head_loc, r1, r2, r3};
}
void testPlot(const std::vector<std::vector<std::vector<double>>> ini) {
    std::ofstream initial_loc_stream("test/no3-initial-loc.txt");
    for(int i = 0; i < N; i++) {
        auto r = ini[i][0];
        initial_loc_stream << r[0] << "\t" << r[1] << "\t" << r[2] << "\n";
    }
    initial_loc_stream.close();
}
void gnuplot() {
    FILE *gp;
	if ((gp = popen(GNUPLOT_PATH, "w")) == NULL) {
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}
    fprintf(gp, "set xlabel \"x\"\n");
    fprintf(gp, "set ylabel \"y\"\n");
    fprintf(gp, "set zlabel \"z\"\n");
    // fprintf(gp, "set lmargin 10\n");
    // fprintf(gp, "set bmargin 4\n");
    // fprintf(gp, "set xlabel font 'Arial,13'\n");
    // fprintf(gp, "set ylabel font 'Arial,13'\n");
    fprintf(gp, "set pointsize 0.5\n");
    // fprintf(gp, "set yrange [0:0.6]\n");
    // fileの読み込みが最後じゃないとうまくいかなかったのなんでだろう。比較的重めのIOなので、それが間接的な原因になってる気がする。コルーチン化したら直る気がするけどめんどいのでこれでいいや。
    fprintf(gp, "sp \"./test/no3-initial-loc.txt\" title \"location\" pt 7\n");
    fflush(gp);
	
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
	fprintf(gp, "exit\n");
	pclose(gp);
}
std::vector<std::vector<double>> normalize(std::vector<std::vector<double>> vs) {
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;
    for(int i = 0; i < N; i++) {
        sum_x += vs[i][0];
        sum_y += vs[i][1];
        sum_z += vs[i][2];
    }
    const double x_mean = sum_x / N;
    const double y_mean = sum_y / N;
    const double z_mean = sum_z / N;
    std::vector<std::vector<double>> normalized = vs;
    for(int i = 0; i < N; i++) {
        normalized[i][0] -= x_mean;
        normalized[i][1] -= y_mean;
        normalized[i][2] -= z_mean;
    }
    return normalized;
}
std::vector<std::vector<double>> randomizedV() {
    std::vector<std::vector<double>> vs(N, std::vector<double>(3));
    for(int i = 0; i < N; i++) {
        vs[i] = vRamdomConfig();
    }
    return normalize(vs);
}
std::vector<std::vector<std::vector<double>>> zip(std::vector<std::vector<double>> rs, std::vector<std::vector<double>> vs) {
    std::vector<std::vector<std::vector<double>>> res(N, std::vector<std::vector<double>>(2));
    if(rs.size() == vs.size() && rs.size() == N) {
        for(int i = 0; i < N; i++) {
            res[i][0] = rs[i];
            res[i][1] = vs[i];
        }
    }
    return res;
}
// ここで初期位置、初期速度、初期温度の設定
std::vector<std::vector<std::vector<double>>> config() {
    std::vector<std::vector<double>> initial_loc(N, std::vector<double>(3));
    // どの分子か->座標か速度か->その値(3次元)
    std::vector<double> head_loc = {0, 0, 0};
    // ここでは正四面体の頭の座標を更新していいく 2*aずつヅラせばいい
    for(int i = 0; i < N / 4; i++) {
        const std::vector<std::vector<double>> hex_loc = culcHexLoc(head_loc);
        initial_loc[i * 4] = hex_loc[0];
        initial_loc[i * 4 + 1] = hex_loc[1];
        initial_loc[i * 4 + 2] = hex_loc[2];
        initial_loc[i * 4 + 3] = hex_loc[3];
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
    std::vector<std::vector<double>> initial_v = randomizedV();
    return zip(initial_loc, initial_v);
}

int main() {
    const std::vector<std::vector<std::vector<double>>> ini = config();
    std::cout << ini.size() << std::endl;
    testPlot(ini);
    gnuplot();
}