#ifndef F_H
#define F_H
#include <cmath>
double func(double x) {
  double y;
  int i;
  y = x;
  for(i = 1;i<=10;i++){
    y = y+sin(x*i)/pow(2,i);
  }
  return y;
}

#endif
