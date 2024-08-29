#include <iostream>
#include <vector>
#include <cmath>

// 计算给定基数 base 的 Halton 序列的第 index 项
double haltonSequence(int index, int base) {
    double result = 0.0;
    double f = 1.0;
    while (index > 0) {
        f = f / base;
        result = result + f * (index % base);
        index = index / base;
    }
    return result;
}

// 生成维度为 dim，长度为 n 的 Halton 序列
std::vector<std::vector<double>> generateHaltonSequence(int n, int dim) {
    std::vector<std::vector<double>> sequence(n, std::vector<double>(dim));
    std::vector<int> bases = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};  // 常用的素数作为基数

    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < dim; ++d) {
            sequence[i][d] = haltonSequence(i + 1, bases[d]);
        }
    }
    
    return sequence;
}

int main() {
    int n = 10;   // 序列长度
    int dim = 3;  // 维度

    std::vector<std::vector<double>> haltonSeq = generateHaltonSequence(n, dim);

    // 输出生成的 Halton 序列
    for (const auto& vec : haltonSeq) {
        for (double val : vec) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
