OUT = ./test/myProgram.out
t = hoge.cpp
comp: ${t}
	g++ ./${t} -o $(OUT) -std=c++11

run: $(OUT)
	$(OUT)