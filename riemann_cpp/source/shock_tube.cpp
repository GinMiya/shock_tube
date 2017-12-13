// 衝撃波管問題(Shock tube test)
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits>
// #include <direct.h>
#include <sys/stat.h>
#include <omp.h>

using namespace std;
#include "riemann_function.hpp"
#define GNUPLOT

void output_ini_cond(quantity *Left, quantity *Right){
  cout << "===========================================" << endl;
  cout << "Initial Condition." << endl;
  cout << "Left(dens, vel, press) = "
  << Left->get_density()  << " "
  << Left->get_velocity() << " "
  << Left->get_pressure() << " "
  << Left->get_sound()    << endl;
  cout << "Left(dens, vel, press) = "
  << Right->get_density()  << " "
  << Right->get_velocity() << " "
  << Right->get_pressure() << " "
  << Right->get_sound()    << endl;
  cout << "===========================================" << endl;
}
void set_initial_condition(quantity *Left, quantity *Right){
  cout << "Left Side:: ";
  Left->set_condition();
  cout << "Right Side:: ";
  Right->set_condition();
  output_ini_cond(Left, Right);
}
void set_gammma_function(double G[8]){
  G[0] = (GAMMA-1.0)/(2.0*GAMMA);
  G[1] = (GAMMA+1.0)/(2.0*GAMMA);
  G[2] = 2.0*GAMMA/(GAMMA-1.0);
  G[3] = 2.0/(GAMMA-1.0);
  G[4] = 2.0/(GAMMA+1.0);
  G[5] = (GAMMA-1.0)/(GAMMA+1.0);
  G[6] = (GAMMA-1.0)/2.0;
  G[7] = GAMMA-1.0;
}

int main(void){
  quantity *W, *WL, *WR, *WS;
  W  = new quantity[1];
  WL = new quantity[1]; WR = new quantity[1]; WS = new quantity[1];
  double gamma_func[8];
  double bin = (double)DOMLEN/CELLS, x_pos = 0.0;
  double u_sound = 0.0;
  double p_star, u_star;
  double t = 0.0;
  int i, nt = 0;
  char filename[50];

  //保存用フォルダの作成
  char dirname[50], _dirname[50];
  sprintf(dirname,"data");
  mkdir(dirname,0777);
  chmod(dirname,0777);

  // 左右の初期物理量設定と確認用出力
  set_initial_condition(WL, WR);
  // γに関する汎用値の計算
  set_gammma_function(gamma_func);
  // 真空状態は扱えないので終了
  double check_vacuum = gamma_func[3]*(WL->get_sound() + WR->get_sound())
                        - (WL->get_velocity() - WR->get_velocity());
  if (check_vacuum < 0.0) {
    cout << "***Vacuum is generates by data***" << endl;
    cout << "***Program stopped***" << endl;
    return 0;
  }

  // 特性線の間(star region)における圧力と速度を求める関数
  STARPU(&p_star, &u_star, WL, WR, WS, gamma_func);
  // 時間発展
  while(t<TIMEOUT){
    strcpy(_dirname,dirname);
		sprintf(filename,"/shock_tube_%d.dat",nt);
    strcat(_dirname,filename);
    ofstream outputfile(_dirname);
		for(i=1; i<=CELLS; i++){
			x_pos = ((double)i-0.5)*bin;
			u_sound = (x_pos-DIAPH)/t;
			SAMPLE(u_sound, WS, W, WL, WR, gamma_func);
      outputfile  << x_pos << "\t"
                  << W->get_density() << "\t"
                  << W->get_velocity() << "\t"
                  << W->get_pressure()/MPA << "\t"
                  << W->get_pressure()/(W->get_density()*gamma_func[7]*MPA) << endl;
		}
    outputfile .close();
		t += dt;
		nt++;
  }

  #ifdef GNUPLOT
    FILE *fpg;
    fpg = popen("gnuplot", "w");
    fputs("load './gnuplot/density_snap.gp'\n", fpg);
    fputs("load './gnuplot/velocity_snap.gp'\n", fpg);
    fputs("load './gnuplot/pressure_snap.gp'\n", fpg);
    fputs("load './gnuplot/energy_snap.gp'\n", fpg);
    fflush(fpg);
    pclose(fpg);
  #endif
  return 0;
}
