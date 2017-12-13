#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "phys_quantity.hpp"
typedef phys_quantity quantity;

//Star Regionにの圧力と速度を求めるための関数
void STARPU(double *P, double *U, quantity *WL, quantity *WR, quantity *WS, double G[8]);
//Star Region内の圧力と速度のおおよその値である圧力PMと速度UMを求める関数//Two rarefaction近似
void GUESSP(quantity *WL, quantity *WR, quantity *WS, double G[8]);
//Riemann Solverでのf関数の計算(左側)
void PREFUNL(double *FL, double *FLD, quantity *WS, quantity *WL, double G[8]);
//Riemann Solverでのf関数の計算(右側)
void PREFUNR(double *FR, double *FRD, quantity *WS, quantity *WR, double G[8]);
//Star Region内での圧力と速度がわかっていて，波のパターンごとの解を求める
void SAMPLE(double S, quantity *WS, quantity *W, quantity *WL, quantity *WR, double G[8]);
