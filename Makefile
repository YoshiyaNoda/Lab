OUT = myProgram.out
T = hoge.cpp
compile: ${T}
	g++ ./${T} -o $(OUT) -std=c++11

test: $(OUT)
	./$(OUT)