/*
A monotone cubic interpolation algorithm 
Author:  Rick Su 
Date:    2024.8.28
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

struct SplineSegment{
    double a, b, c, d, x;// a + b(x-x_i) + c(x-x_i)^2 + d(x-x_i)^3
};

//compute the gradient of every part
std::vector<double> computeGradient(const std::vector<double>& x, const std::vector<double>& y){
    size_t n = x.size();
    std::vector<double> gradient(n-1);
    for(size_t i = 0; i < n-1; i++){
        gradient[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
    }
    return gradient;
}

//compute the second derivative of every part,ensure the monotonicity
std::vector<double> computeSecondDerivative(const std::vector<double>& x, const std::vector<double>& y){
    size_t n = x.size();
    std::vector<double> gradient = computeGradient(x, y);
    std::vector<double> second_derivative(n);
    //solve the first point
    second_derivative[0] = 0;
    for(size_t i = 1; i < n-1; i++){
        double s = (gradient[i-1] * (x[i+1] - x[i]) + gradient[i] * (x[i] - x[i-1])) / (x[i+1] - x[i-1]);
        second_derivative[i] = std::max(std::min(2.0, s), -2.0);
    }
    //solve the last point
    second_derivative[n-1] = gradient[n-2];
    return second_derivative;
}

//build the monotone cubic spline
std::vector<SplineSegment> buildMonotoneCubicSpline(const std::vector<double>& x, const std::vector<double>& y){
    size_t n = x.size();
    std::vector<SplineSegment> spline(n-1);
    std::vector<double> gradient = computeGradient(x, y);
    std::vector<double> second_derivative = computeSecondDerivative(x, y);
    for(size_t i = 0; i < n-1; i++){
        double h = x[i+1] - x[i];
        spline[i].a = y[i];
        spline[i].b = second_derivative[i];
        spline[i].c = (3.0 * (y[i+1]-y[i])/h-2.0*second_derivative[i]-second_derivative[i+1])/h;
        spline[i].d = (2.0 * (y[i]-y[i+1]) + second_derivative[i] + second_derivative[i+1]) / (h*h);
    }
    return spline;
}

//compute the value of the monotone cubic spline
double interpolate(const std::vector<SplineSegment>& spline, double x){
   if (spline.empty()) {
    throw std::invalid_argument("spline is empty");
   }

   // Find the segment containing x
   auto it = std::lower_bound(spline.begin(), spline.end(), x, [](const SplineSegment& s, double x) { return s.x < x; });
   if (it == spline.end()) {
    it = spline.begin();
   }else if (it == spline.end()){
    it = spline.end();
   }else{
    --it;
   }

   double dx =  x - it->x;
   return it->a + (it->b + (it->c + it->d * dx) * dx) * dx;
}

int main(){
    //inputs
    std::vector<double> x = {0, 1, 2, 3, 4, 5};
    std::vector<double> y = {0, 1, 0, 1, 0, 1};
    //build the monotone cubic spline
    std::vector<SplineSegment> spline = buildMonotoneCubicSpline(x, y);
    //interpolate
    double x_interp=2.5;
    double y_interp=interpolate(spline, x_interp);
    std::cout<<"The interpolated value at x="<<x_interp<<" is "<<y_interp<<std::endl;
    return 0;
}