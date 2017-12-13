#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//状況設定(キーボード入力の場合のみ使用text読み込みの場合はいらない)
#define DIAPH	0.5			//初期接触不連続面の位置
#define DOMLEN	1.0			//Tubeの長さ
#define CELLS	1000		//binの個数(解像度)
#define GAMMA	1.4			//比熱比
#define MPA		1.0			//規格化定数
#define TIMEOUT	0.01		//出力時間
#define dt		0.00078125	//タイムステップ

//関数定義
//Star Region内の圧力と速度のおおよその値である圧力PMと速度UMを求める関数//Two rarefaction近似
void GUESSP(double WL[4], double WR[4], double WS[4], double G[8]){
	//変数定義
	double PPV;			//pressure based on primitive variables(4.47)式参照
	double PMIN,PMAX;	//左右の圧力を比較し，大小を比較するためのもの
	double QMAX;		//圧力比
	double QUSER = 2.0;	//スイッチングパラメータ
	double PQ,PTL,PTR;	//two rarefaction approximationの際に使うもの
	double GEL,GER;		//two shock appriximationの際に使うもの
	
	//compute guess pressure from PVRS Riemann solver
	PPV = ((WL[2]+WR[2])/2.0)+((WL[1]-WR[1])*(WL[0]+WR[0])*(WL[3]+WR[3])/8.0);	//PPVの計算(4.47)式
	PPV = fmax(0.0, PPV);		//0.0とPPVの大きいほうをとる
	PMIN = fmin(WL[2], WR[2]);	//左右の小さいほうをとる
	PMAX = fmax(WL[2], WR[2]);	//左右の大きいほうをとる
	QMAX = PMAX/PMIN;			//圧力比
	
	if(QMAX<=QUSER && (PMIN<=PPV && PPV<=PMAX)){
		//select PVRS Riemann solver
		WS[2] = PPV;		//Star Regionでの圧力はPPVとする
	}
	else if(PPV<PMIN){
		//Two Rarefanction approximation (4.46)式
		PQ = pow((WL[2]/WR[2]),G[0]);
		WS[1] = (PQ*WL[1]/WL[3]+WR[1]/WR[3]+G[3]*(PQ-1.0))/(PQ/WL[3]+1.0/WR[3]);
		PTL = 1.0+G[6]*(WL[1]-WS[1])/WL[3];
		PTR = 1.0+G[6]*(WS[1]-WR[1])/WR[3];
		WS[2] =0.5*(pow(WL[2]*PTL,G[2])+pow(WR[2]*PTR,G[2]));
		printf("Two-Rarefaction approximation selected.\n");
//		WS[2] = pow(((WL[3]+WR[3]-G[6]*(WR[1]-WL[1]))/(pow(WL[3]/WL[2],G[0])+pow(WR[3]/WR[2],G[0]))),G[2]);
	}
	else{
		//select two-shock Riemann solver with PVRS as estimate (4.48)式 これが採用されるはず
		GEL = sqrt((G[4]/WL[0])/(G[5]*WL[2]+PPV));
		GER = sqrt((G[4]/WR[0])/(G[5]*WR[2]+PPV));
		WS[2] = (GEL*WL[2]+GER*WR[2]-(WR[1]-WL[1]))/(GEL+GER);
		printf("Two-Shock approximation selected.\n");
	}
}

//Riemann Solverでのf関数の計算(左側)
void PREFUNL(double *FL, double *FLD, double WS[4], double WL[4], double G[8]){
	double AL,BL,PRATL,QRTL;	//定数A_K,B_KとFのためのもの
	if (WS[2]<WL[2]){
		//rarefaction wave
		PRATL = WS[2]/WL[2];
		*FL = G[3]*WL[3]*(pow(PRATL,G[0]) - 1.0);
		*FLD = (1.0/(WL[0]*WL[3]))/pow(PRATL,G[1]);
	}else{
		//shock wave
		AL = G[4]/WL[0];
		BL = G[5]*WL[2];
		QRTL = sqrt(AL/(BL + WS[2]));
		*FL = (WS[2] - WL[2])*QRTL;
		*FLD = (1.0 - (WS[2] - WL[2])/(2.0*(BL + WS[2])))*QRTL;
	}
}
//Riemann Solverでのf関数の計算(右側)
void PREFUNR(double *FR, double *FRD, double WS[4], double WR[4], double G[8]){
	double AR,BR,PRATR,QRTR;
	if (WS[2]<WR[2]){
		//rarefaction wave
		PRATR = WS[2]/WR[2];
		*FR= G[3]*WR[3]*(pow(PRATR,G[0]) - 1.0);
		*FRD = (1.0/(WR[0]*WR[3]))/pow(PRATR,G[1]);
	}else{
		//shock wave
		AR = G[4]/WR[0];
		BR = G[5]*WR[2];
		QRTR = sqrt(AR/(BR + WS[2]));
		*FR = (WS[2] - WR[2])*QRTR;
		*FRD = (1.0 - (WS[2] - WR[2])/(2.0*(BR + WS[2])))*QRTR;
	}
}

//Star Regionにの圧力と速度を求めるための関数
void STARPU(double *P, double *U, double WL[4], double WR[4], double WS[4], double G[8]){
	int i;
	int NRITER = 10;					//反復回数
	double TOLPRE = 1/pow(10.0,6.0);	//p0であるという条件
	double UDIFF;						//左右波の速度差
	double FL,FLD,FR,FRD;				//左右の圧力に関する関数
	double CHANGE;						//圧力の差分
	
	//Guessed value PSTART is computed
	GUESSP(WL, WR ,WS, G);					//反復回数が減るような圧力の推察
	printf("guess result PM=%lf\n",WS[2]);	//推定圧力の出力
	
	UDIFF = WR[1] - WL[1];		//波の左右での速度差計算
	
	printf("===========================================\n");
	printf("Interation number:Change\n");
	printf("===========================================\n");
	for(i=1; i<=NRITER; i++)
	{
		PREFUNL(&FL, &FLD, WS, WL, G);				//関数Fの計算(左側)
		PREFUNR(&FR, &FRD, WS, WR, G);				//関数Fの計算(右側)
		*P = WS[2]-(FL+FR+UDIFF)/(FLD+FRD);			//新たな圧力の計算(Newton法)
		CHANGE = 2.0*fabs((*P-WS[2]))/(*P+WS[2]);	//圧力の変化量を計算
		printf("i=%d, CHANGE=%lf\n",i,CHANGE);		//何回目のinterationかと圧力変化量の出力
		if(CHANGE<TOLPRE){
			continue;								//変化量が10^-6よりも小さいときはinterationしない
		}
		if(*P<0.0){
			*P=TOLPRE;								//圧力が負になったらそのときは最も小さい10^-6とする
		}
		WS[2]=(*P);
	}												//interationの終了
	
		printf("Divergence in Newton-Raphson interation\n");
		*U = (WL[1] + WR[1] + FR - FL)/2.0;			//interationでわかったstar regionでの圧力を使って速度を計算
		WS[1]=(*U);
		//Compute velocity in star region
		printf("===========================================\n");
		printf("Pressure\t Velocity\n");
		printf("===========================================\n");
		printf("%.10f\t%.10f\n",*P/MPA,*U);
		printf("===========================================\n\n");
}

//Star Region内での圧力と速度がわかっていて，波のパターンごとの解を求める
void SAMPLE(double S, double WS[4], double W[4], double WL[4], double WR[4], double G[4]){
	double SHL,STL;		//Left rarefaction head,tail
	double SHR,STR;		//Right rarefaction head,tail
	double CML,CMR;		//rarefaction内の音速の計算
	double PML,PMR;		//SL,SRを計算するためのもの
	double SL,SR;		//左右のshockの速度
	
	if (S<WS[1]){
//		printf("check:first brench->");
		//Sampling point lies to the left of the contact discontinuity
		if(WS[2]<WL[2]){
//			printf("check:left wave->");
			//left rarefaction
			SHL=WL[1]-WL[3];	//左側rarefactionの速度
			if(S<SHL){
//				printf("check:left rarefaction in head.\n");
				//samples point is lecft data state
				W[0] = WL[0];
				W[1] = WL[1];
				W[2] = WL[2];
			}else{
				CML = WL[3]*pow((WS[2]/WL[2]),G[0]);
				STL = WS[1] - CML;
				if(S>STL){
//					printf("check:left rarefaction in tail.\n");
					//samples point is star left state
					W[0] = WL[0]*pow((WS[2]/WL[2]),(1.0/GAMMA));
					W[1] = WS[1];
					W[2] = WS[2];
				}else{
//					printf("check:left fan in rarefaction\n");
					//sampled point is inside left fan
					W[1] = G[4]*(WL[3] + G[6]*WL[1] + S);
					W[3] = G[4]*(WL[3] + G[6]*(WL[1] - S));
					W[0] = WL[0]*pow((W[3]/WL[3]),G[3]);
					W[2] = WL[2]*pow((W[3]/WL[3]),G[2]);
				}
			}
		}else{
			//left shock
//			printf("check:left shock.\n");
			PML = WS[2]/WL[2];
			SL = WL[1] - WL[3]*sqrt(G[1]*PML + G[0]);
			if(S<SL){
				//sampled point is left data state
//				printf("check:left shock out.\n");
				W[0] = WL[0];
				W[1] = WL[1];
				W[2] = WL[2];
			}else{
				//sampled point is star left state
//				printf("check:left shock in.\n");
				W[0] = WL[0]*(PML + G[5])/(PML*G[5] + 1.0);
				W[1] = WS[1];
				W[2] = WS[2];
			}
		}
	}else if(WS[2]>WR[2]){
//		printf("check:right wave->");
		PMR = WS[2]/WR[2];
		SR = WR[1] + WR[3]*sqrt(G[1]*PMR + G[0]);
		if(S>SR){
			//Sampled point is right data state
//			printf("check:right shock out.\n");
			W[0] = WR[0];
			W[1] = WR[1];
			W[2] = WR[2];
		}else{
			//Sampled point is Star Right state
//			printf("check:right shock in.\n");
			W[0] = WR[0]*(PMR + G[5])/(PMR*G[5] + 1.0);
			W[1] = WS[1];
			W[2] = WS[2];
		}
	}else{
		//right rarefaction
//		printf("check:right rarefaction->");
		SHR = WR[1]+WR[3];	//head wave
		if(S>SHR){
			//sampled point is right data state
//			printf("check:right rarefaction in head.\n");
			W[0] = WR[0];
			W[1] = WR[1];
			W[2] = WR[2];
		}else{
			CMR = WR[3]*pow((WS[2]/WR[2]),G[0]);
			STR = WS[1] + CMR;
			if(S<STR){
				//sampled point is star right state
//				printf("check:right rarefaction in tail.\n");
				W[0] = WR[0]*pow((WS[2]/WR[2]),(1.0/GAMMA));
				W[1] = WS[1];
				W[2] = WS[2];
			}else{
			//Sampled point is inside right fan
//				printf("check:right fan in rarefaction.\n");
				W[1] = G[4]*(-WR[3] + G[6]*WR[1] + S);
				W[3] = G[4]*(WR[3] - G[6]*(WR[1] - S));
				W[0] = WR[0]*pow((W[3]/WR[3]),G[3]);
				W[2] = WR[2]*pow((W[3]/WR[3]),G[2]);
			}
		}
	}
}

int main (void){
	//時間計測の開始
	clock_t start,end;
	start=clock();
	printf("starttime[ms]:%d\n",start);
	
	//変数定義
	double W[4],WL[4],WR[4],WS[4];	//物理量の配列(密度[0]，速度[1]，圧力[2]，音速[3])
	double G[8];			//γの値とγによる関数の定義
	double S;					//衝撃波の速度
	double DX;					//binの長さ
	double XPOS;				//x座標
	double P,U;					//Star Regionでの圧力，速度を求めるためのポインタ扱い
	double t;					//時間
	int i,t2=0;					//ループ用変数
	
	char filename[50];
	
	//配列の初期化
	for(i=0;i<4;i++){
		W[i]=0.0;
		WL[i]=0.0;
		WR[i]=0.0;
		WS[i]=0.0;
	}
//	printf("check:初期化\n");
	
	//Set Up条件をinitial_data.datから読み込む
/*	double DIAPH;		//初期接触不連続面の位置
	double DOMLEN;		//Tubeの長さ
	double CELLS;		//binの個数(解像度)
	double GAMMA;		//比熱比
	double TIMEOUT;		//出力時間
	double MPA;			//規格化定数
	FILE *fp;
	fp= fopen("initial_data.dat","r");
	//初期データ読み込み
	if((fp == NULL))
	{
		printf("初期データを読み込めませんでした．\n");
		return 1;
	}
	fscanf(fp,"%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&DOMLEN,&DIAPH,&CELLS,&GAMMA,&TIMEOUT,&(WL[0]),&(WL[1]),&(WL[2]),&(WR[0]),&(WR[1]),&(WR[2]),&MPA);
		//変数確認用標準出力
		printf("=================================\n");
		printf("Initial Condition.\n");
		printf("DOMLEN=%lf\nDIAPH=%lf\nCELLS=%d\nGAMMA=%lf\nTIMEOUT=%d\nDL=%lf\nUL=%lf\nPL=%lf\nDR=%lf\nUR=%lf\nPR=%lf\nMPA=%lf\n",DOMLEN,DIAPH,CELLS,GAMMA,TIMEOUT,WL[0],WL[1],WL[2],WR[0],WR[1],WR[2],MPA);
		printf("=================================\n");
	fclose(fp);
*/	
	//左側密度，速度，圧力の入力
	printf("Left Side 'density' 'velocity' 'pressure'=");
	scanf("%lf%lf%lf",&(WL[0]),&(WL[1]),&(WL[2]));
	//右側密度，速度，圧力の入力
	printf("Right Side 'density' 'velocity' 'pressure'=");
	scanf("%lf%lf%lf",&(WR[0]),&(WR[1]),&(WR[2]));
	
	//変数確認用標準出力
	printf("===========================================\n");
	printf("Initial Condition.\n");
	printf("(DL,UL,PL)=(%lf,%lf,%lf)\n",WL[0],WL[1],WL[2]);
	printf("(DR,UR,PR)=(%lf,%lf,%lf)\n",WR[0],WR[1],WR[2]);
	printf("DOMLEN=%lf\nDIAPH=%lf\nCELLS=%d\nGAMMA=%lf\nTIMEOUT=%lf\n",DOMLEN,DIAPH,CELLS,GAMMA,TIMEOUT);
	printf("===========================================\n\n\n");
	
	
	//GAMMAの関数の計算
	G[0] = (GAMMA-1.0)/(2.0*GAMMA);
	G[1] = (GAMMA+1.0)/(2.0*GAMMA);
	G[2] = 2.0*GAMMA/(GAMMA-1.0);
	G[3] = 2.0/(GAMMA-1.0);
	G[4] = 2.0/(GAMMA+1.0);
	G[5] = (GAMMA-1.0)/(GAMMA+1.0);
	G[6] = (GAMMA-1.0)/2.0;
	G[7] = GAMMA-1.0; 
	
	//音速の計算
	WL[3] = sqrt(GAMMA*WL[2]/WL[0]);
	WR[3] = sqrt(GAMMA*WR[2]/WR[0]);
	
	//真空は扱わないのでプログラム終了
	if (G[3]*(WL[3]+WR[3]) < (WR[1]-WL[1]))
	{
		//the initial data is such that vacuum is generated.program stopped.
		printf("***Vacuum is generates by data***\n***Program stopped***\n");
		return 1;
	}
	
	//Star Regionでの圧力と速度を求める関数の実行
	STARPU(&P,&U,WL,WR,WS,G);
	
	//bin幅の計算
	DX = DOMLEN/(double)CELLS;
	

	//TIMEOUT時の解を求める
	while(t<=TIMEOUT){
		FILE *fp;
		sprintf(filename,"exact_%d.dat",t2);
		fp = fopen(filename,"w");
		for(i=1; i<=CELLS; i++){
			XPOS = ((double)i-0.5)*DX;	//x座標の計算(何個目のbinか)
			S = (XPOS-DIAPH)/t;			//そのxでのTIMEOUT時における音速の計算
			
			//solution at point (x,t)=(XPOS-DIAPH,TIMEPUT) is found
			SAMPLE(S,WS,W,WL,WR,G);
			//exact solution prifiles are written to exact.dat(x座標と密度，速度，圧力，エネルギーの出力)
			fprintf(fp, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", XPOS, W[0], W[1], W[2]/MPA, W[2]/(W[0]*G[7]*MPA));
//			printf("FILE OUT:count %d\tXPOS=%lf\n",i,XPOS);
		}
		t = t + dt;
		t2++;
		fclose(fp);
	}
	
	/*計測時間の終了*/
	end=clock();
	printf("endtime:%d[ms]\n",end);
	printf("uptime:%d[ms]\n",end-start);
	
	return 0;
	
}
