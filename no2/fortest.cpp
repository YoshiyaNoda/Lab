#include<iostream>
#include<vector>
std::vector<double> divideBy6(std::vector<double> v) {
    std::vector<double> res = v;
    for(int i = 0; i < res.size(); i++) {
        res[i] = res[i] / 6; // 多分これ6.0と書かないとダメ
    }
    return res;
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
int main() {
    const std::vector<double> hoge1 = {6,12,6};
    const std::vector<double> hoge2 = {6,12,6};
    const std::vector<double> hoge3 = {6,12,6};
    const std::vector<double> hoge = divideBy6(vectorsSum({hoge1, hoge2, hoge3}));
    for(int i = 0; i < hoge.size(); i++) {
        std::cout << hoge[i] << std::endl;
    }
}