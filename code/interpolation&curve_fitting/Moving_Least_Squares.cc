#include <iostream>
#include <vector>
#include <cmath>

// Gaussian kernel function to calculate weights
double gaussianKernel(double distance, double bandwidth) {
    return exp(- (distance * distance) / (2.0 * bandwidth * bandwidth));
}

// Moving Least Squares fitting function
double movingLeastSquares(const std::vector<double>& x, const std::vector<double>& y, double queryPoint, int degree, double bandwidth) {
    int n = x.size();
    
    // Construct weight matrix
    std::vector<double> weights(n);
    for (int i = 0; i < n; ++i) {
        weights[i] = gaussianKernel(x[i] - queryPoint, bandwidth);
    }

    // Construct X matrix and Y vector
    std::vector<std::vector<double>> X(n, std::vector<double>(degree + 1));
    std::vector<double> WY(n);
    for (int i = 0; i < n; ++i) {
        double xi = 1.0;
        for (int j = 0; j <= degree; ++j) {
            X[i][j] = xi * weights[i];
            xi *= x[i] - queryPoint;  // Polynomial basis function
        }
        WY[i] = y[i] * weights[i];
    }

    // Calculate (X^T * W * X) matrix
    std::vector<std::vector<double>> XT_W_X(degree + 1, std::vector<double>(degree + 1, 0.0));
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k < n; ++k) {
                XT_W_X[i][j] += X[k][i] * X[k][j];
            }
        }
    }

    // Calculate (X^T * W * Y) vector
    std::vector<double> XT_W_Y(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        for (int k = 0; k < n; ++k) {
            XT_W_Y[i] += X[k][i] * WY[k];
        }
    }

    // Solve the linear system (XT_W_X) * a = XT_W_Y using simple Gaussian elimination
    std::vector<double> a(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        int maxRow = i;
        for (int k = i + 1; k <= degree; ++k) {
            if (fabs(XT_W_X[k][i]) > fabs(XT_W_X[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(XT_W_X[i], XT_W_X[maxRow]);
        std::swap(XT_W_Y[i], XT_W_Y[maxRow]);

        for (int k = i + 1; k <= degree; ++k) {
            double factor = XT_W_X[k][i] / XT_W_X[i][i];
            for (int j = i; j <= degree; ++j) {
                XT_W_X[k][j] -= factor * XT_W_X[i][j];
            }
            XT_W_Y[k] -= factor * XT_W_Y[i];
        }
    }

    for (int i = degree; i >= 0; --i) {
        a[i] = XT_W_Y[i] / XT_W_X[i][i];
        for (int k = i - 1; k >= 0; --k) {
            XT_W_Y[k] -= XT_W_X[k][i] * a[i];
        }
    }

    // Return the fitted result, i.e., the value at queryPoint
    double result = 0.0;
    double xi = 1.0;
    for (int j = 0; j <= degree; ++j) {
        result += a[j] * xi;
        xi *= queryPoint;
    }

    return result;
}

int main() {
    // Input data points
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {1.0, 2.0, 0.0, 2.0, 1.0};

    // Set parameters
    double queryPoint = 2.5;
    int degree = 2;         // Degree of the polynomial to fit
    double bandwidth = 1.0; // Bandwidth, controls the smoothness

    // Calculate the fitted value at the queryPoint
    double result = movingLeastSquares(x, y, queryPoint, degree, bandwidth);

    std::cout << "Fitted value at x = " << queryPoint << " is y = " << result << std::endl;

    return 0;
}
