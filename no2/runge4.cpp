#include <iostream>
#include <fstream>
#include <cmath>
#define GNUPLOT_PATH "gnuplot"
#include <unistd.h>


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

auto f_A = [](double arg_x_a, double arg_x_b, double arg_x_c) -> double {
    return k / m * (arg_x_b - 2 * arg_x_a);
};
auto f_B = [](double arg_x_a, double arg_x_b, double arg_x_c) -> double {
    return k / m * (arg_x_a + arg_x_c - 2 * arg_x_b);
};
auto f_C = [](double arg_x_a, double arg_x_b, double arg_x_c) -> double {
    return k / m * (arg_x_b + L * w - 2 * arg_x_c);
};
double nextHalfX(double v, double x) {
    return x + (dt / 2) * v;
}
double nextHalfV(double v, double x_a, double x_b, double x_c, std::function<double(double, double, double)> f) {
    return v + (dt / 2) * f(x_a, x_b, x_c);
}
double nextV(double v, double x_a, double x_b, double x_c, std::function<double(double, double, double)> f) {
    return v + dt * f(x_a, x_b, x_c);
}
double nextX(double x, double v) {
    return x + v * dt;
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
    double half_v_a;
    double half_v_b;
    double half_v_c;
    double half_x_a;
    double half_x_b;
    double half_x_c;
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
        half_x_a = nextHalfX(v_a, x_a);
        half_x_b = nextHalfX(v_b, x_b);
        half_x_c = nextHalfX(v_c, x_c);
        half_v_a = nextHalfV(v_a, x_a, x_b, x_c, f_A);
        half_v_b = nextHalfV(v_b, x_a, x_b, x_c, f_B);
        half_v_c = nextHalfV(v_c, x_a, x_b, x_c, f_C);
        x_a = nextX(x_a, half_v_a);
        x_b = nextX(x_b, half_v_b);
        x_c = nextX(x_c, half_v_c);
        v_a = nextV(v_a, half_x_a, half_x_b, half_x_c, f_A);
        v_b = nextV(v_b, half_x_a, half_x_b, half_x_c, f_B);
        v_c = nextV(v_c, half_x_a, half_x_b, half_x_c, f_C);
    }
    potential.close();
    kinetic.close();
    total.close();
    gnuplot();
    return 0;
}