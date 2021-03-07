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
double nextHalfX(double v, double x) {
    return x + (dt / 2) * v;
}
double nextHalfV(double v, std::vector<double> xs, std::function<double(std::vector<double>)> f) {
    return v + (dt / 2) * f(xs);
}
double nextV(double v, std::vector<double> vs, std::vector<double> xs, std::function<double(std::vector<double>)> f) {
    const std::vector<double> half_xs = {nextHalfX(vs[0], xs[0]), nextHalfX(vs[1], xs[1]), nextHalfX(vs[2], xs[2])};
    return v + dt * f(half_xs);
}
double nextX(double v, double x, std::vector<double> xs, std::function<double(std::vector<double>)> f) {
    const double half_v = nextHalfV(v, xs, f);
    return x + half_v * dt;
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
    // fileの読み込みが最後じゃないとうまくいかなかったのなんでだろう。結構重めのIOなので、それが間接的な原因になってる気がする。コルーチン化したら直る気がするけどめんどいのでこれでいいや。
    fprintf(gp, "plot \"./test/potential.txt\" title \"potential\" pt 7, \"./test/kinetic.txt\" title \"kinetic\" pt 9, \"./test/total_energy.txt\" title \"total\" pt 11\n");
    fflush(gp);
	
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
	fprintf(gp, "exit\n");
	pclose(gp);
}

double culcU(double x_a, double x_b, double x_c) {
    return 0.5 * k * ((x_a - L) * (x_a - L)) + 0.5 * k * ((x_b - x_a - L) * (x_b - x_a - L)) + 0.5 * k * ((x_c - x_b - L) * (x_c - x_b - L)) + 0.5 * k * ((x_c - L * (w - 1)) * (x_c - L * (w - 1)));
}
double culcK(double v_a, double v_b, double v_c) {
    return 0.5 * m * (v_a * v_a + v_b * v_b + v_c * v_c);
}

int main() {
    std::ofstream potential("test/potential.txt");
    std::ofstream kinetic("test/kinetic.txt");
    std::ofstream total("test/total_energy.txt");

    double v_a = v_a0;
    double v_b = v_b0;
    double v_c = v_c0;
    double x_a = x_a0;
    double x_b = x_b0;
    double x_c = x_c0;
    double t = t_min;

    int N = (t_max - t_min) / dt + 1; // iもintなので同じ型にしておく
    for(int i = 0; i < N; i++) {
        const double U = culcU(x_a, x_b, x_c);
        const double K = culcK(v_a, v_b, v_c);
        potential << t << "\t" << U << "\n";
        kinetic << t << "\t" << K << "\n";
        total << t << "\t" << K + U << "\n";

        // 値の更新
        t += dt;
        const std::vector<double> xs = {x_a, x_b, x_c};
        const std::vector<double> vs = {v_a, v_b, v_c};
        x_a = nextX(v_a, x_a, xs, f_A);
        x_b = nextX(v_b, x_b, xs, f_B);
        x_c = nextX(v_c, x_c, xs, f_C);
        v_a = nextV(v_a, vs, xs, f_A);
        v_b = nextV(v_b, vs, xs, f_B);
        v_c = nextV(v_c, vs, xs, f_C);
    }
    potential.close();
    kinetic.close();
    total.close();
    gnuplot();
    return 0;
}