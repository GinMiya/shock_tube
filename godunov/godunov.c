#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//#define GAMMA	1.4		//比熱比
#define GAMMA	(5.0/3.0)	//比熱比
#define DOMLEN	1.0		//tubeの長さ
#define DIAPH	0.5		//仕切りの位置
#define CELLS	100		//セルの数
#define TEND	0.002	//終了時間
#define MPA		1.0		//規格化定数
#define CFLCOE	0.9		//CFL係数

//関数定義
//定数出力
void PRINT(){
	//値の標準出力
	printf("=======================================\n");
	printf("Dates of This Simulation\n");
	printf("=======================================\n");
	printf("CFLCOE=%lf\n",CFLCOE);
	printf("DOMLEN=%lf\n",DOMLEN);
	printf("DIAPH=%lf\n",DIAPH);
	printf("CELLS=%d\n",CELLS);
	printf("GAMMA=%lf\n",GAMMA);
}
//GAMMAに関する関数
void GFUNC(double G[8]){
	G[0] = (GAMMA-1.0)/(2.0*GAMMA);
	G[1] = (GAMMA+1.0)/(2.0*GAMMA);
	G[2] = 2.0*GAMMA/(GAMMA-1.0);
	G[3] = 2.0/(GAMMA-1.0);
	G[4] = 2.0/(GAMMA+1.0);
	G[5] = (GAMMA-1.0)/(GAMMA+1.0);
	G[6] = (GAMMA-1.0)/2.0;
	G[7] = GAMMA-1.0;
}
//真空の条件判定
int VACUUM(double INIWL[4], double INIWR[4], double G[8]){
	//真空は扱わないのでプログラム終了
	if(G[3]*(INIWL[3]+INIWR[3])<(INIWR[1]-INIWL[1])){
		//the initial data is such that vacuum is generated.program stopped.
		printf("***Vacuum is generates by data***\n***Program stopped***\n");
		return 1;
	}//end if
}//end VACUUM
//初期物理量の出力
void INICOND(double INIWL[4], double INIWR[4], double G[8]){
	//音速の計算と出力⇒CFL条件に使う
	INIWL[3] = sqrt(GAMMA*INIWL[2]/INIWL[0]);
	INIWR[3] = sqrt(GAMMA*INIWR[2]/INIWR[0]);
	printf("Initial Condition\n");
	printf("(DL,UL,PL,AL)=(%lf,%lf,%lf,%lf)\n",INIWL[0],INIWL[1],INIWL[2],INIWL[3]);
	printf("(DR,UR,PR,AR)=(%lf,%lf,%lf,%lf)\n",INIWR[0],INIWR[1],INIWR[2],INIWR[3]);
	printf("=============================================================================\n\n");
}
//初期状態の設定関数CELLに物理量を入れるのと初期化/
void INITIAL(double U[][3], double INIWL[4], double INIWR[4], double DX, double G[8]){
	int i,j,k;
	//保存量の初期化
	for(i=0;i<CELLS+2;i++){
		for(k=0;k<3;k++){
			U[i][k]=0.0;
		}//end for k
	}//end for i
	//境界抜きの物理量を入れる(保存形)/初期は半分のセルを境界にして物理量の分配(衝撃波管問題のため)
	for(i=1;i<DIAPH/DX+1;i++){
		//密度
		U[i][0]=INIWL[0];
		//密度×速度
		U[i][1]=INIWL[0]*INIWL[1];
		//単位質量あたりの全エネルギー
		U[i][2]=INIWL[2]/G[7]+0.5*INIWL[0]*INIWL[1]*INIWL[1];
	}//end for i
	for(i=DIAPH/DX+1;i<CELLS+1;i++){
		//密度
		U[i][0]=INIWR[0];
		//密度×速度
		U[i][1]=INIWR[0]*INIWR[1];
		//単位質量あたりの全エネルギー
		U[i][2]=INIWR[2]/G[7]+0.5*INIWR[0]*INIWR[1]*INIWR[1];
	}//end for i
}//end INITIAL
//境界条件関数
void BOUNDI(double U[][3]){
	int k;
	for(k=0;k<3;k++){
		U[0][k]=0.0;
		U[CELLS+1][k]=0.0;
		U[0][k]=U[1][k];
		U[CELLS+1][k]=U[CELLS][k];
	}//end for
}//end BOUNDI
//CFL条件
void CFLCOND(double U[][3], double *DT, double DX, double G[8]){
	int i,j,k;
	double SMAX=-pow(10.0,6.0),SN;
	for(i=1;i<CELLS+1;i++){
		SN=sqrt(GAMMA*fabs(U[i][2]-0.5*U[i][0]*U[i][1]*U[i][1]/(U[i][0]*U[i][0]))*G[7]/U[i][0])+(U[i][1]/U[i][0]);
		if(fabs(SN)>=SMAX){
			SMAX=fabs(SN);
		}//end if
	}//end for
	*DT=CFLCOE*DX/SMAX;
}//end CFLCOND
/********************************
*ここから下はRiemann問題である *
********************************/
void GUESSP(double WL[4], double WR[4], double WS[4], double G[8]){
	//変数定義
	double PPV;			//pressure based on primitive variables(4.47)式参照
	double PMIN,PMAX;	//左右の圧力を比較し，大小を比較するためのもの
	double QMAX;		//圧力比
	double QUSER = 2.0;	//スイッチングパラメータ
	double PQ,PTL,PTR;	//two rarefaction approximationの際に使うもの
	double GEL,GER;		//two shock appriximationの際に使うもの
	//compute guess pressure from PVRS Riemann solver
	PPV = (0.5*(WL[2]+WR[2])/2.0)+(0.125*(WL[1]-WR[1])*(WL[0]+WR[0])*(WL[3]+WR[3]));	//PPVの計算(4.47)式
	PPV = fmax(0.0, PPV);		//0.0とPPVの大きいほうをとる
	PMIN = fmin(WL[2], WR[2]);	//左右の小さいほうをとる
	PMAX = fmax(WL[2], WR[2]);	//左右の大きいほうをとる
	QMAX = PMAX/PMIN;			//圧力比
	//select two-shock Riemann solver with PVRS as estimate (4.48)式 これが採用されるはず
		GEL = sqrt((G[4]/WL[0])/(G[5]*WL[2]+PPV));
		GER = sqrt((G[4]/WR[0])/(G[5]*WR[2]+PPV));
		WS[2] = (GEL*WL[2]+GER*WR[2]-(WR[1]-WL[1]))/(GEL+GER);
}//end GUESSP
//Riemann Solverでのf関数の計算(左側)
void PREFUNL(double *FL, double *FLD, double WS[4], double WL[4], double G[8]){
	double AL,BL,PRATL,QRTL;	//定数A_K,B_KとFのためのもの
	if (WS[2]<=WL[2]){
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
}//end OREFUNL
//Riemann Solverでのf関数の計算(右側)
void PREFUNR(double *FR, double *FRD, double WS[4], double WR[4], double G[8]){
	double AR,BR,PRATR,QRTR;
	if (WS[2]<=WR[2]){
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
}//end PREFUNR
//Star Regionでの圧力と速度を求めるための関数
void STARPU(double WL[4], double WR[4], double WS[4], double G[8]){
	int j;
	int NRITER = 10;					//反復回数
	double TOLPRE = 1.0/pow(10.0,6.0);	//p0であるという条件
	double UDIFF;						//左右波の速度差
	double FL,FLD,FR,FRD;				//左右の圧力に関する関数
	double CHANGE;						//圧力の差分
	double PS,US;						//StarRegionの圧力と速度の推定値
	//Guessed value PSTART is computed
	GUESSP(WL,WR,WS,G);				//反復回数が減るような圧力の推察
	UDIFF = WR[1] - WL[1];			//波の左右での速度差計算
	for(j=1; j<=NRITER; j++){
		PREFUNL(&FL, &FLD, WS, WL, G);				//関数Fの計算(左側)
		PREFUNR(&FR, &FRD, WS, WR, G);				//関数Fの計算(右側)
		PS = WS[2]-(FL+FR+UDIFF)/(FLD+FRD);			//新たな圧力の計算(Newton法)
		CHANGE = 2.0*fabs((PS-WS[2]))/(PS+WS[2]);	//圧力の変化量を計算
		if(CHANGE<TOLPRE){
			continue;		//変化量が10^-6よりも小さいときはinterationしない
		}
		if(PS<0.0){
			PS=TOLPRE;		//圧力が負になったらそのときは最も小さい10^-6とする
		}
		WS[2]=PS;
	}						//interationの終了
	//interationでわかったstar regionでの圧力を使って速度を計算
	US = (WL[1] + WR[1] + FR - FL)/2.0;
	WS[1]=US;
}//end STARPU
//局所Riemann問題
void RIEMANN(double WS[4], double W[4], double WL[4], double WR[4], double G[8], double XPOS){
	double SHL,STL;		//Left rarefaction head,tail
	double SHR,STR;		//Right rarefaction head,tail
	double CML,CMR;		//rarefaction内の音速の計算
	double PML,PMR;		//SL,SRを計算するためのもの
	double SL,SR;		//左右のshockの速度
	//printf("RIEMANN pressure velocity in StarRegion=%lf\t%lf\n",WS[2],WS[1]);
	if(WS[1]>0.0){
		if(WS[2]>WL[2]){
			//shock wave
			PML=WS[2]/WL[2];
			//左のshock速度
			SL = WL[1]-WL[3]*sqrt(G[1]*PML+G[0]);
			if(SL>0.0){
				//patern(a1):Outside of Shock 				printf("patern(a1)\n");
				W[0]=WL[0];
				W[1]=WL[1];
				W[2]=WL[2];
				W[3]=WL[3];
			}else if(SL<=0.0){
				//patern(a2):StarRegion 				printf("patern(a2)\n");
				W[0]=WL[0]*pow((WS[2]/WL[2]),(1.0/(GAMMA)));
				W[1]=WS[1];
				W[2]=WS[2];
				W[3]=WL[3]*pow((WS[2]/WL[2]),(1.0/G[2]));
			}
		}else{
			//rarefaction wave
			SHL=WL[1]-WL[3];	//左側rarefactionの速度
			CML=WL[3]*pow((WS[2]/WL[2]),G[0]);//a*Lの計算(4.55)
			STL=WS[1]-CML;	//左側rarefactionの速度
			if((SHL>0.0)&&(STL>0.0)){
				//patern(a3):Outside of Rarefaction 				printf("patern(a3)\n");
				W[0]=WL[0];
				W[1]=WL[1];
				W[2]=WL[2];
				W[3]=WL[3];
			}else if((SHL<0.0)&&(STL<0.0)){
				//patern(a4):StarRegion 				printf("patern(a4)\n");
				W[0]=WL[0]*pow((WS[2]/WL[2]),(1.0/(GAMMA)));
				W[1]=WS[1];
				W[2]=WS[2];
				W[3]=WL[3]*pow((WS[2]/WL[2]),(1.0/G[2]));
			}else{
				//patern(a5):In fan of rarefaction 				printf("patern(a5)\n");
				W[1]=G[4]*(WL[3]+G[6]*WL[1]);
				W[3]=G[4]*(WL[3]+G[6]*WL[1]);//=W[1]
				W[0]=WL[0]*pow((W[3]/WL[3]),G[3]);
				W[2]=WL[2]*pow((W[3]/WL[3]),G[2]);
			}//end patern(a5)
		}//end rarefaction
	}else if(WS[1]<=0.0){
		if(WS[2]>WR[2]){
			//shock wave
			PMR=WS[2]/WR[2];
			SR=WR[1]+WR[3]*sqrt(G[1]*PMR+G[0]);
			if(SR<0.0){
				//patern(b1) printf("patern(b1) at %lf\n",XPOS);
				W[0]=WR[0];
				W[1]=WR[1];
				W[2]=WR[2];
				W[3]=WR[3];
			}else if(SR>=0.0){
				//patern(b2) 				printf("patern(b2) at %lf\n",XPOS);
				W[0]=WR[0]*((PMR+G[5])/(G[5]*PMR+1.0));
				W[1]=WS[1];
				W[2]=WS[2];
				W[3]=WR[3]*pow((WS[2]/WR[2]),G[0]);
			}
		}else{
			//rarefaction wave
			SHR=WR[1]+WR[3];	//右側rarefaction headの速度
			CMR=WR[3]*pow((WS[2]/WR[2]),G[0]);	//a*Rの計算(4.62)
			STR=WS[1]+CMR;		//右側rarefaction tailの速度
			if((SHR<0.0)&&(STR<0.0)){
				//patern(b3) 				printf("patern(b3) at %lf\n",XPOS);
				W[0]=WR[0];
				W[1]=WR[1];
				W[2]=WR[2];
				W[3]=WR[2];
			}else if((SHR>0.0)&&(STR>0.0)){
				//patern(b4) 				printf("patern(b4) at %lf\n",XPOS);
				W[0]=WR[0]*pow((WS[2]/WR[2]),(1.0/GAMMA));
				W[1]=WS[1];
				W[2]=WS[2];
				W[3]=CMR;
			}else{
				//patern(b5) 				printf("patern(b5) at %lf\n",XPOS);
				W[1]=G[4]*(-WR[3]+G[6]*WR[1]);
				W[3]=G[4]*(WR[3]-G[6]*WR[1]);//=-W[1]
				W[0]=WR[0]*pow((W[3]/WR[3]),G[3]);
				W[2]=WR[2]*pow((W[3]/WR[3]),G[2]);
			}//end patern(a5)
		}//end rarefaction
	}else{
		//printf("!!error no patern!! at %lf\n",XPOS);
	}//end else
}//end RIEMANN
//フラックス計算
void FLUXES(double U[][3], double FLUX[][3], double WL[4], double WR[4], double W[4], double WS[4], double G[8], double DX, double XPOS){
	int i,j,k;
	for(i=0;i<CELLS+1;i++){
		//Step4-1:FLUXの初期化
		for(k=0;k<3;k++){
			W[k]=0.0;
			WS[k]=0.0;
			FLUX[i][k]=0.0;
		}//end for
		//Step4-2:座標ごとのRieamann問題を解いてFluxを求める
		XPOS=(i+0.5)*DX;
		//保存量を速度，圧力について書きなおしてやる必要がある
			//Riemann問題用左側密度,速度,圧力,音速
			WL[0]=U[i][0];
			WL[1]=U[i][1]/(U[i][0]);
			WL[2]=(U[i][2]-0.5*WL[0]*WL[1]*WL[1])*G[7];
			WL[3]=sqrt(fabs(GAMMA*WL[2]/WL[0]));
			//Riemann問題用右側密度,速度,圧力,音速
			WR[0]=U[i+1][0];
			WR[1]=U[i+1][1]/(U[i+1][0]);
			WR[2]=(U[i+1][2]-0.5*WR[0]*WR[1]*WR[1])*G[7];
			WR[3]=sqrt(fabs(GAMMA*WR[2]/WR[0]));
		//Riemann問題を解くにあたってのStarRegionでの圧力,速度の見積もり
		STARPU(WL,WR,WS,G);
		//RIEMANNの実行で密度,速度,圧力をもとめる
		RIEMANN(WS,W,WL,WR,G,XPOS);
		//Eular方程式のフラックス
		FLUX[i][0]=W[0]*W[1];
		FLUX[i][1]=W[0]*W[1]*W[1]+W[2];
		FLUX[i][2]=W[1]*0.5*(G[2]*W[2]+W[0]*W[1]*W[1]);
	}//end for
}//end FLUX
//次のタイムステップへアップデート
void UPDATE(double DTODX, double U[][3], double FLUX[][3]){
	int i,j,k;
	for(i=1;i<CELLS+1;i++){
		U[i][0] = U[i][0]+DTODX*(FLUX[i-1][0]-FLUX[i][0]);
		U[i][1] = U[i][1]+DTODX*(FLUX[i-1][1]-FLUX[i][1]);
		U[i][2] = U[i][2]+DTODX*(FLUX[i-1][2]-FLUX[i][2]);
	}//end for
}//end UPDATE

//main関数
int main(void){
	//定数定義
	int i,j,k,N=1;
	double G[8];						//GAMMA関数
	double INIWL[4],INIWR[4];			//初期物理量配列
	double U[CELLS+2][3],FLUX[CELLS+2][3];	//保存量とフラックス
	double W[4],WS[4];					//Riemann問題の解とStarRegionの物理量
	double WL[4],WR[4];					//Riemann問題の左右の物理量
	double FL,FLD,FR,FRD;				//Riemann問題用の圧力関数
	double r,u,p,e;						//出力用密度,速度,圧力,内部エネルギー
	double P_,U_;						//STARPUのためのポインタ
	double DX=0.0,DT=0.0,DTODX=0.0;		//セル幅,タイムステップ,DT/DX
	double XPOS;						//x軸の位置
	double TIME=0.0;					//時間
	//ファイル出力準備
	char dataname[50];
	FILE *fp;FILE *fp0;
	fp0=fopen("godunov_0.dat","w");
	//GAMMAに関する関数の定義
	GFUNC(G);
	//各定数の出力
	PRINT();
	//初期物理量の入力
		printf("=============================================================================\n");
		//左側密度，速度，圧力の入力
		printf("Left  Side 'density' 'velocity' 'pressure'=");
		scanf("%lf%lf%lf",&(INIWL[0]),&(INIWL[1]),&(INIWL[2]));
		//右側密度，速度，圧力の入力
		printf("Right Side 'density' 'velocity' 'pressure'=");
		scanf("%lf%lf%lf",&(INIWR[0]),&(INIWR[1]),&(INIWR[2]));
		printf("=============================================================================\n");
	//真空判定
	VACUUM(INIWL,INIWR,G);
	//初期物理量の確認
	INICOND(INIWL,INIWR,G);
	//セルの長さ/初期保存量の格納
	DX=DOMLEN/((double)CELLS);
	INITIAL(U,INIWL,INIWR,DX,G);
	//初期物理量の出力
	for(i=1;i<CELLS+1;i++){
		//保存量を物理量で書く
		r=U[i][0];
		u=U[i][1]/r;
		p=(U[i][2]-0.5*r*u*u)*G[7];
		e=p/(r*G[7]*MPA);
		fprintf(fp0,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",FLUX[i][0],FLUX[i][1],FLUX[i][2],U[i][0],U[i][1],U[i][2]);
		//fprintf(fp0,"%lf\t%lf\t%lf\t%lf\t%lf\n",i*DX,r,u,p,e);
	}//end for
	fclose(fp0);
	//時間ループ
	do{
		//境界条件
		BOUNDI(U);
		////CFL条件からタイムステップ決定
		CFLCOND(U,&DT,DX,G);
		DTODX=DT/DX;
		////時間をDTだけ進める
		TIME+=DT;
		////物理量を計算,出力
		printf("=============================================================================\n");
		printf("TIME=%lf\tNo.%d TIMESTEP DT=%lf\t DT/DX=%lf\n",TIME,N,DT,DTODX);
		printf("-----------------------------------------------------------------------------\n");
		////フラックスの計算
		FLUXES(U,FLUX,WL,WR,W,WS,G,DX,XPOS);
		////保存量のアップデート
		UPDATE(DTODX,U,FLUX);
		sprintf(dataname,"godunov_%lf.dat",TIME);
		fp=fopen(dataname,"w");
		for(i=1;i<CELLS+1;i++){
			//保存量を物理量で書く
			r=U[i][0];
			u=U[i][1]/r;
			p=(U[i][2]-0.5*r*u*u)*G[7];
			e=p/(r*G[7]*MPA);
			//fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",FLUX[i][0],FLUX[i][1],FLUX[i][2],U[i][0],U[i][1],U[i][2]);
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",(i-0.5)*DX,r,u,p,e);
		}//end for
		fclose(fp);
		N++;
/*		if(N>200){
			return 0;
		}//end if*/
	}while(TIME<TEND);//end while
	return 0;
}//end main
