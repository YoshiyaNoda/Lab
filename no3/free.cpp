// https://www.fml.t.u-tokyo.ac.jp/~izumi/CMS/MD/models.pdf を参考にした
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#define GNUPLOT_PATH "gnuplot"
#include <vector>

constexpr int N = 256;
constexpr int M = 3; //step数
constexpr double dt = 0.0050;
constexpr double m = 1.0;
constexpr double sigma = 1.0;
constexpr double yps = 1.0;
// const double a = std::pow(4.0 / 1.6, 1.0 / 3.0);
const double a = std::pow(4.0 / 0.1, 1.0 / 3.0);
constexpr int rand_min = 100;
constexpr int rand_max = -100;

std::random_device rd;
std::default_random_engine eng(rd());
std::uniform_real_distribution<double> distr(rand_min, rand_max);
auto norm3d = [](std::vector<double> v) -> double {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
};

std::vector<double> vRamdomConfig() {
    return {distr(eng), distr(eng), distr(eng)};
}
std::vector<std::vector<double>> culcHexLoc(std::vector<double> head_loc) {
    const std::vector<double> r1 = {head_loc[0] + a, head_loc[1] + a, head_loc[2]};
    const std::vector<double> r2 = {head_loc[0], head_loc[1] + a, head_loc[2] + a};
    const std::vector<double> r3 = {head_loc[0] + a, head_loc[1], head_loc[2] + a};
    return {head_loc, r1, r2, r3};
}
void testPlot(std::vector<std::vector<std::vector<double>>> ini) {
    std::ofstream initial_loc_stream("test/no3-initial-loc.txt");
    for(int i = 0; i < N; i++) {
        auto r = ini[i][0];
        initial_loc_stream << r[0] << "\t" << r[1] << "\t" << r[2] << "\n";
    }
    initial_loc_stream.close();
}
void testVPlot(std::vector<std::vector<double>> vs_data) {
    std::ofstream initial_v_stream("test/no3-free-velocity.txt");
    for(int i = 0; i < vs_data.size(); i++) {
        const auto vs = vs_data[i];
        auto v = std::sqrt(vs[0] * vs[0] + vs[1] * vs[1] + vs[2] * vs[2]);
        initial_v_stream << i << "\t" << v << "\n";
    }
    initial_v_stream.close();
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
    fprintf(gp, "set pointsize 0.5\n");
    fprintf(gp, "sp \"./test/no3-initial-loc.txt\" title \"location\" pt 7\n");
    fflush(gp);
	
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
	fprintf(gp, "exit\n");
	pclose(gp);
}
void gnuplotV() {
    FILE *gp;
	if ((gp = popen(GNUPLOT_PATH, "w")) == NULL) {
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}
    fprintf(gp, "set xlabel \"t\"\n");
    fprintf(gp, "set ylabel \"v\"\n");
    fprintf(gp, "set pointsize 0.5\n");
    fprintf(gp, "plot \"./test/no3-free-velocity.txt\" title \"velocity\" pt 7\n");
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
double updateK(std::vector<std::vector<std::vector<double>>> data) {
    double k = 0.0;
    for(int i = 0; i < N; i++) {
        const std::vector<double> v = data[i][1];
        k += 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    return k;
}
double updateU(std::vector<std::vector<std::vector<double>>> data) {
    double p = 0.0;
    for(int i = 0; i < N; i++) {
        const std::vector<double> r_i = data[i][0];
       for(int j = 0; j < i; j++) {
           const std::vector<double> r_j = data[j][0];
           const double sq_r = (r_j[0] - r_i[0]) * (r_j[0] - r_i[0]) + (r_j[1] - r_i[1]) * (r_j[1] - r_i[1]) + (r_j[2] - r_i[2]) * (r_j[2] - r_i[2]);
           p += 4 * yps * (std::pow((sigma * sigma) / sq_r, 6) - std::pow((sigma * sigma) / sq_r, 3));
       }
    }

    return p;
}
std::vector<double> culcF(std::vector<double> r_i, std::vector<std::vector<std::vector<double>>> erased) {
    std::vector<double> f = {0.0, 0.0, 0.0};
    const int size = erased.size();
    for(int j = 0; j < size; j++) {
        const std::vector<double> r_j = erased[j][0];
        const double sq_r = (r_j[0] - r_i[0]) * (r_j[0] - r_i[0]) + (r_j[1] - r_i[1]) * (r_j[1] - r_i[1]) + (r_j[2] - r_i[2]) * (r_j[2] - r_i[2]);
        f[0] += 24 * yps * (2 * std::pow((sigma * sigma) / sq_r, 6) - std::pow((sigma * sigma) / sq_r, 3)) / sq_r * (r_j[0] - r_i[0]);
        f[1] += 24 * yps * (2 * std::pow((sigma * sigma) / sq_r, 6) - std::pow((sigma * sigma) / sq_r, 3)) / sq_r * (r_j[1] - r_i[1]);
        f[2] += 24 * yps * (2 * std::pow((sigma * sigma) / sq_r, 6) - std::pow((sigma * sigma) / sq_r, 3)) / sq_r * (r_j[2] - r_i[2]);
    }
    return f;
}
std::vector<double> culcH(std::vector<double> v_i, std::vector<double> f_i) {
    std::vector<double> h = {0.0, 0.0, 0.0};
    h[0] = v_i[0] + (1.0 / (2.0 * m)) * dt * f_i[0];
    h[1] = v_i[1] + (1.0 / (2.0 * m)) * dt * f_i[1];
    h[2] = v_i[2] + (1.0 / (2.0 * m)) * dt * f_i[2];
    return h;
}
std::vector<std::vector<std::vector<double>>> updateMData(std::vector<std::vector<std::vector<double>>> data) {
    std::vector<std::vector<std::vector<double>>> res = data;
    std::vector<double> h = {0.0, 0.0, 0.0}; // 部分速度
    for(int i = 0; i < N; i++) {
        const std::vector<double> r_i = data[i][0];
        const std::vector<double> v_i = data[i][1];
        std::vector<std::vector<std::vector<double>>> erased = data;
        erased.erase(erased.begin() + i);
        const std::vector<double> f_i = culcF(r_i, erased);
        // if(i == 1) {
        //     std::cout << f_i[0] << "\t" << f_i[1] << "\t" << f_i[2] << std::endl;
        // }
        h = culcH(v_i, f_i);
        // std::cout << "h\t" << norm3d(h) << std::endl;
        // 座標更新
        res[i][0] = {r_i[0] + dt * h[0], r_i[1] + dt * h[1], r_i[2] + dt * h[2]};
    }
    for(int i = 0; i < N; i++) {
        std::vector<std::vector<std::vector<double>>> erased = data;
        erased.erase(erased.begin() + i);
        const std::vector<double> r_i = data[i][0];
        const std::vector<double> v_i = data[i][1];
        const std::vector<double> f_i = culcF(r_i, erased);
        h = culcH(v_i, f_i);
        const std::vector<double> r_i_next = res[i][0];
        std::vector<std::vector<std::vector<double>>> erased_next = res;
        erased_next.erase(erased_next.begin() + i);
        const std::vector<double> f_i_next = culcF(r_i_next, erased_next);
        // 速度更新
        res[i][1] = {h[0] + (1.0 / (2.0 * m)) * dt * f_i_next[0], h[1] + (1.0 / (2.0 * m)) * dt * f_i_next[1], h[2] + (1.0 / (2.0 * m)) * dt * f_i_next[2]};
    }
    return res;
}
int main() {
    std::vector<std::vector<std::vector<double>>> m_data = config();
    // testPlot(ini);
    std::ofstream kinetic("test/no3-free-kinetic.txt");
    std::ofstream potential("test/no3-free-potential.txt");
    std::ofstream energy("test/no3-free-energy.txt");
    std::ofstream velocity("test/no3-free-velocity.txt");
    
    double k = 0.0;
    double p = 0.0;
    double t = 0.0;
    std::vector<std::vector<double>> vs_data(M, std::vector<double>(3));
    for(int i = 0; i < M; i++) {
        k = updateK(m_data);
        p = updateU(m_data);
        auto old = m_data; //
        m_data = updateMData(m_data);
        // for(int j = 0; j < m_data.size(); j++) {
        //     auto before = old[j][1];
        //     auto after = m_data[j][1];

        //     if(norm3d(before) > norm3d(after)) {
        //         std::cout << "index!" << "\t" << j << "\t" << "  before: " << norm3d(before) << "  after: " << norm3d(after) << std::endl;
        //     }
        // }
        // vs_data[i] = m_data[1][1]; //速度変化観察用
        std::cout << i << "\t" << k + p << "\t" << p << "\t" << k << std::endl;
        kinetic << t << "\t" << k << "\n";
        potential << t << "\t" << p << "\n";
        energy << t << "\t" << k + p << "\n";
        // velocity << t << "\t" << m_data[2][1][0] << "\t" << m_data[2][1][1] << "\t" << m_data[2][1][2] << "\t" << "\n";
        t += dt;
    }
    // testVPlot(vs_data);//速度テスト用
    // gnuplotV(); 
    testPlot(m_data);
    gnuplot();
    kinetic.close();
    potential.close();
    energy.close();
    // velocity.close();
}