#include <iostream> 
#define GNUPLOT_PATH "gnuplot"

int main(void) {
    FILE *gp;	// For gnuplot
	// gnuplotの起動コマンド
	if ((gp = popen(GNUPLOT_PATH, "w")) == NULL) {	// gnuplotをパイプで起動
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}

	// --- gnuplotにコマンドを送る --- //
	fprintf(gp, "set xrange [-10:10]\n"); // 範囲の指定(省略可)
	fprintf(gp, "set yrange [-1:1]\n");

	fprintf(gp, "plot sin(x)\n"); 	//sin(x)を描く
	
	fflush(gp); // バッファに格納されているデータを吐き出す（必須）
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
	fprintf(gp, "exit\n"); // gnuplotの終了
	pclose(gp);

    return 0;
}