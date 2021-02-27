#include <iostream>

using namespace std;

int main(void) {
    FILE *gp;
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "plot sin(x)");
    pclose(gp);
    return 0;
}