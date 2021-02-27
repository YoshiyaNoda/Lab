#include <iostream>
using namespace std;

// 定数
const float m = 0.7;
const float k = 0.2;
const float g = 9.81;
const float dt = 0.1;

// 初期条件
const float t0 = 0;
const float v0 = 0;

float nextV(float v) {
    return (g - k/m * v) * dt + v;
}

int main() {
    float v = v0;
    for(int i = 0; i < 20; i++) {
        v = nextV(v);
        cout << v << endl;
    }
    return 0;
}
