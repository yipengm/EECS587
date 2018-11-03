#ifndef F_H
#define F_H

#include <cmath>

double f(double x) {
    double ans = 0;
    for (int i = 100; i >= 1; --i) {
        double temp = 0;
        for (int j = i; j >= 1; --j) temp += pow(x + 0.5*j, -3.3);
        ans += sin(x + temp) / pow(1.3, i);
    }
    return ans;
}


#endif