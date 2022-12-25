#include "SlowSolution.h"
#include "OptimizedSolution.h"

int main() {
    int num = 100;
    int Nmax = 100000;
    double eps = 0.00000001;
    std::cout << "Slow solution\n";
    slowSolution(num, num, Nmax, eps);
    std::cout << "Optimized solution:\n";
    optimizedSolution(num, num, Nmax, eps);
    return 0;
}
