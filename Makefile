OUT = ./test/myProgram.out
T = hoge.cpp
comp: ${T}
	g++ ./${T} -o $(OUT) -std=c++11

run: $(OUT)
	$(OUT)