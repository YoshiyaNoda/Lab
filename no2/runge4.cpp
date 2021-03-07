#include <iostream>
#include <fstream>
#include <cmath>
#define GNUPLOT_PATH "gnuplot"
#include <vector>


// 定数
constexpr double m = 1.0;
constexpr double k = 1.0;
constexpr double L = 1.0;
constexpr double w = 5.0;

// 初期条件
constexpr double v_a0 = 0;
constexpr double v_b0 = 0;
constexpr double v_c0 = 0;
constexpr double x_a0 = 1.0;
constexpr double x_b0 = 2.0;
constexpr double x_c0 = 3.0;

// その他条件
constexpr double dt = 0.01;
constexpr double t_min = 0;
constexpr double t_max = 100;

auto f_A = [](std::vector<double> xs) -> double {
    return k / m * (xs[1] - 2 * xs[0]);
};
auto f_B = [](std::vector<double> xs) -> double {
    return k / m * (xs[0] + xs[2] - 2 * xs[1]);
};
auto f_C = [](std::vector<double> xs) -> double {
    return k / m * (xs[1] + L * w - 2 * xs[2]);
};
double tempX(double v, double x, double stride) {
    return x + stride * v;
}
double tempV(double v, double stride, std::vector<double> xs, std::function<double(std::vector<double>)> f) {
    return v + stride * f(xs);
}
std::vector<double> tempXs(std::vector<double> vs, std::vector<double> xs, double stride) {
    return {tempX(vs[0], xs[0], stride), tempX(vs[1], xs[1], stride), tempX(vs[2], xs[2], stride)};
}
std::vector<double> tempVs(std::vector<double> vs, std::vector<double> xs, double stride) {
    return {tempV(vs[0], stride, xs, f_A), tempV(vs[1], stride, xs, f_B), tempV(vs[2], stride, xs, f_C)};
}
std::vector<double> vectorSum(std::vector<double> a, std::vector<double> b) {
    //とりあえずめんどくさいので長さは同じ前提
    std::vector<double> res = a;
    for(int i = 0; i < a.size(); i++) {
        res[i] += b[i];
    }
    return res;
}
std::vector<double> vectorsSum(std::vector<std::vector<double>> vectors) {
    std::vector<double> res(vectors[0].size());
    for(int i = 0; i < vectors.size(); i++) {
        res = vectorSum(res, vectors[i]);
    }
    return res;
}
std::vector<double> devideBy6(std::vector<double> v) {
    std::vector<double> res = v;
    for(int i = 0; i < res.size(); i++) {
        res[i] = res[i] / 6;
    }
    return res;
}
std::vector<std::vector<double>> update(std::vector<std::vector<double>> d) {
    const std::vector<double> vs_before = d[0];
    const std::vector<double> xs_before = d[1];

    const std::vector<double> vs_temp_1 = tempVs(vs_before, xs_before, dt / 2);
    const std::vector<double> xs_temp_1 = tempXs(vs_before, xs_before, dt / 2);

    const std::vector<double> vs_temp_2 = tempVs(vs_temp_1, xs_temp_1, dt / 2);
    const std::vector<double> xs_temp_2 = tempXs(vs_temp_1, xs_temp_1, dt / 2);

    const std::vector<double> vs_temp_3 = tempVs(vs_temp_2, xs_temp_2, dt);
    const std::vector<double> xs_temp_3 = tempXs(vs_temp_2, xs_temp_2, dt);

    const std::vector<double> afterVs = devideBy6(vectorsSum({tempVs(vs_before, xs_before, dt), tempVs(vs_temp_1, xs_temp_1, 2 * dt), tempVs(vs_temp_2, xs_temp_2, 2 * dt), tempVs(vs_temp_3, xs_temp_3, dt)}));
    const std::vector<double> afterXs = devideBy6(vectorsSum({tempXs(vs_before, xs_before, dt), tempXs(vs_temp_1, xs_temp_1, 2 * dt), tempXs(vs_temp_2, xs_temp_2, 2 * dt), tempXs(vs_temp_3, xs_temp_3, dt)}));
    return {afterVs, afterXs};
}
double culcU(std::vector<double> xs) {
    const double x_a = xs[0];
    const double x_b = xs[1];
    const double x_c = xs[2];
    return 0.5 * k * ((x_a - L) * (x_a - L)) + 0.5 * k * ((x_b - x_a - L) * (x_b - x_a - L)) + 0.5 * k * ((x_c - x_b - L) * (x_c - x_b - L)) + 0.5 * k * ((x_c - L * (w - 1)) * (x_c - L * (w - 1)));
}
double culcK(std::vector<double> vs) {
    const double v_a = vs[0];
    const double v_b = vs[1];
    const double v_c = vs[2];
    return 0.5 * m * (v_a * v_a + v_b * v_b + v_c * v_c);
}

void gnuplot() {
    FILE *gp;
	if ((gp = popen(GNUPLOT_PATH, "w")) == NULL) {
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}
    fprintf(gp, "set xlabel \"time [s]\"\n");
    fprintf(gp, "set ylabel \"energy [J]\"\n");
    fprintf(gp, "set lmargin 10\n");
    fprintf(gp, "set bmargin 4\n");
    fprintf(gp, "set xlabel font 'Arial,13'\n");
    fprintf(gp, "set ylabel font 'Arial,13'\n");
    fprintf(gp, "set pointsize 0.5\n");
    // fileの読み込みが最後じゃないとうまくいかなかったのなんでだろう。比較的重めのIOなので、それが間接的な原因になってる気がする。コルーチン化したら直る気がするけどめんどいのでこれでいいや。
    fprintf(gp, "plot \"./test/potential.txt\" title \"potential\" pt 7, \"./test/kinetic.txt\" title \"kinetic\" pt 9, \"./test/total_energy.txt\" title \"total\" pt 11\n");
    fflush(gp);
	
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
	fprintf(gp, "exit\n");
	pclose(gp);
}

int main() {
    std::ofstream potential("test/potential.txt");
    std::ofstream kinetic("test/kinetic.txt");
    std::ofstream total("test/total_energy.txt");

    const std::vector<double> xs = {x_a0, x_b0, x_c0};
    const std::vector<double> vs = {v_a0, v_b0, v_c0};
    std::vector<std::vector<double>> d = {vs, xs};
    double t = t_min;

    int N = (t_max - t_min) / dt + 1; // iもintなので同じ型にしておく
    for(int i = 0; i < N; i++) {
        const double U = culcU(d[1]);
        const double K = culcK(d[0]);
        potential << t << "\t" << U << "\n";
        kinetic << t << "\t" << K << "\n";
        total << t << "\t" << K + U << "\n";

        // 値の更新
        t += dt;
        d = update(d);
    }
    potential.close();
    kinetic.close();
    total.close();
    gnuplot();
    return 0;
}