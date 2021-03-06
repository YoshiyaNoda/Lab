#include <iostream>
#include <fstream>
#include <cmath>
#define GNUPLOT_PATH "gnuplot"

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

double nextV_a(double v_a, double x_a, double x_b) {
    return v_a - k / m * dt * (2 * x_a - x_b);
}
double nextV_b(double v_b, double x_a, double x_b, double x_c) {
    return v_b - k / m * dt * (2 * x_b - x_a - x_c);
}
double nextV_c(double v_c, double x_b, double x_c) {
    return v_c - k / m * dt * (2 * x_c - x_b - L * w);
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

    fprintf(gp, "plot 'potential.txt' title 'potential' pt 7, 'kinetic.txt' title 'kinetic' pt 9, 'total_energy.txt' title 'total' pt 11\n");
	
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
    std::ofstream potential("potential.txt");
    std::ofstream kinetic("kinetic.txt");
    std::ofstream total("total_energy.txt");

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
        std::cout << t << "\t" << x_a << std::endl; 
        potential << t << "\t" << U << "\n";
        kinetic << t << "\t" << K << "\n";
        total << t << "\t" << K + U << "\n";

        // 値の更新
        t += dt;
        v_a = nextV_a(v_a, x_a, x_b);
        v_b = nextV_b(v_b, x_a, x_b, x_c);
        v_c = nextV_c(v_c, x_b, x_c);
        x_a = nextX(x_a, v_a);
        x_b = nextX(x_b, v_b);
        x_c = nextX(x_c, v_c);
    }
    potential.close();
    kinetic.close();
    total.close();
    gnuplot();
    return 0;
}