#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//�󋵐ݒ�(�L�[�{�[�h���͂̏ꍇ�̂ݎg�ptext�ǂݍ��݂̏ꍇ�͂���Ȃ�)
#define DIAPH	0.5			//�����ڐG�s�A���ʂ̈ʒu
#define DOMLEN	1.0			//Tube�̒���
#define CELLS	1000		//bin�̌�(�𑜓x)
#define GAMMA	1.4			//��M��
#define MPA		1.0			//�K�i���萔
#define TIMEOUT	0.01		//�o�͎���
#define dt		0.00078125	//�^�C���X�e�b�v

//�֐���`
//Star Region���̈��͂Ƒ��x�̂����悻�̒l�ł��鈳��PM�Ƒ��xUM�����߂�֐�//Two rarefaction�ߎ�
void GUESSP(double WL[4], double WR[4], double WS[4], double G[8]){
	//�ϐ���`
	double PPV;			//pressure based on primitive variables(4.47)���Q��
	double PMIN,PMAX;	//���E�̈��͂��r���C�召���r���邽�߂̂���
	double QMAX;		//���͔�
	double QUSER = 2.0;	//�X�C�b�`���O�p�����[�^
	double PQ,PTL,PTR;	//two rarefaction approximation�̍ۂɎg������
	double GEL,GER;		//two shock appriximation�̍ۂɎg������
	
	//compute guess pressure from PVRS Riemann solver
	PPV = ((WL[2]+WR[2])/2.0)+((WL[1]-WR[1])*(WL[0]+WR[0])*(WL[3]+WR[3])/8.0);	//PPV�̌v�Z(4.47)��
	PPV = fmax(0.0, PPV);		//0.0��PPV�̑傫���ق����Ƃ�
	PMIN = fmin(WL[2], WR[2]);	//���E�̏������ق����Ƃ�
	PMAX = fmax(WL[2], WR[2]);	//���E�̑傫���ق����Ƃ�
	QMAX = PMAX/PMIN;			//���͔�
	
	if(QMAX<=QUSER && (PMIN<=PPV && PPV<=PMAX)){
		//select PVRS Riemann solver
		WS[2] = PPV;		//Star Region�ł̈��͂�PPV�Ƃ���
	}
	else if(PPV<PMIN){
		//Two Rarefanction approximation (4.46)��
		PQ = pow((WL[2]/WR[2]),G[0]);
		WS[1] = (PQ*WL[1]/WL[3]+WR[1]/WR[3]+G[3]*(PQ-1.0))/(PQ/WL[3]+1.0/WR[3]);
		PTL = 1.0+G[6]*(WL[1]-WS[1])/WL[3];
		PTR = 1.0+G[6]*(WS[1]-WR[1])/WR[3];
		WS[2] =0.5*(pow(WL[2]*PTL,G[2])+pow(WR[2]*PTR,G[2]));
		printf("Two-Rarefaction approximation selected.\n");
//		WS[2] = pow(((WL[3]+WR[3]-G[6]*(WR[1]-WL[1]))/(pow(WL[3]/WL[2],G[0])+pow(WR[3]/WR[2],G[0]))),G[2]);
	}
	else{
		//select two-shock Riemann solver with PVRS as estimate (4.48)�� ���ꂪ�̗p�����͂�
		GEL = sqrt((G[4]/WL[0])/(G[5]*WL[2]+PPV));
		GER = sqrt((G[4]/WR[0])/(G[5]*WR[2]+PPV));
		WS[2] = (GEL*WL[2]+GER*WR[2]-(WR[1]-WL[1]))/(GEL+GER);
		printf("Two-Shock approximation selected.\n");
	}
}

//Riemann Solver�ł�f�֐��̌v�Z(����)
void PREFUNL(double *FL, double *FLD, double WS[4], double WL[4], double G[8]){
	double AL,BL,PRATL,QRTL;	//�萔A_K,B_K��F�̂��߂̂���
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
//Riemann Solver�ł�f�֐��̌v�Z(�E��)
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

//Star Region�ɂ̈��͂Ƒ��x�����߂邽�߂̊֐�
void STARPU(double *P, double *U, double WL[4], double WR[4], double WS[4], double G[8]){
	int i;
	int NRITER = 10;					//������
	double TOLPRE = 1/pow(10.0,6.0);	//p0�ł���Ƃ�������
	double UDIFF;						//���E�g�̑��x��
	double FL,FLD,FR,FRD;				//���E�̈��͂Ɋւ���֐�
	double CHANGE;						//���͂̍���
	
	//Guessed value PSTART is computed
	GUESSP(WL, WR ,WS, G);					//�����񐔂�����悤�Ȉ��͂̐��@
	printf("guess result PM=%lf\n",WS[2]);	//���舳�͂̏o��
	
	UDIFF = WR[1] - WL[1];		//�g�̍��E�ł̑��x���v�Z
	
	printf("===========================================\n");
	printf("Interation number:Change\n");
	printf("===========================================\n");
	for(i=1; i<=NRITER; i++)
	{
		PREFUNL(&FL, &FLD, WS, WL, G);				//�֐�F�̌v�Z(����)
		PREFUNR(&FR, &FRD, WS, WR, G);				//�֐�F�̌v�Z(�E��)
		*P = WS[2]-(FL+FR+UDIFF)/(FLD+FRD);			//�V���Ȉ��͂̌v�Z(Newton�@)
		CHANGE = 2.0*fabs((*P-WS[2]))/(*P+WS[2]);	//���͂̕ω��ʂ��v�Z
		printf("i=%d, CHANGE=%lf\n",i,CHANGE);		//����ڂ�interation���ƈ��͕ω��ʂ̏o��
		if(CHANGE<TOLPRE){
			continue;								//�ω��ʂ�10^-6�����������Ƃ���interation���Ȃ�
		}
		if(*P<0.0){
			*P=TOLPRE;								//���͂����ɂȂ����炻�̂Ƃ��͍ł�������10^-6�Ƃ���
		}
		WS[2]=(*P);
	}												//interation�̏I��
	
		printf("Divergence in Newton-Raphson interation\n");
		*U = (WL[1] + WR[1] + FR - FL)/2.0;			//interation�ł킩����star region�ł̈��͂��g���đ��x���v�Z
		WS[1]=(*U);
		//Compute velocity in star region
		printf("===========================================\n");
		printf("Pressure\t Velocity\n");
		printf("===========================================\n");
		printf("%.10f\t%.10f\n",*P/MPA,*U);
		printf("===========================================\n\n");
}

//Star Region���ł̈��͂Ƒ��x���킩���Ă��āC�g�̃p�^�[�����Ƃ̉������߂�
void SAMPLE(double S, double WS[4], double W[4], double WL[4], double WR[4], double G[4]){
	double SHL,STL;		//Left rarefaction head,tail
	double SHR,STR;		//Right rarefaction head,tail
	double CML,CMR;		//rarefaction���̉����̌v�Z
	double PML,PMR;		//SL,SR���v�Z���邽�߂̂���
	double SL,SR;		//���E��shock�̑��x
	
	if (S<WS[1]){
//		printf("check:first brench->");
		//Sampling point lies to the left of the contact discontinuity
		if(WS[2]<WL[2]){
//			printf("check:left wave->");
			//left rarefaction
			SHL=WL[1]-WL[3];	//����rarefaction�̑��x
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
	//���Ԍv���̊J�n
	clock_t start,end;
	start=clock();
	printf("starttime[ms]:%d\n",start);
	
	//�ϐ���`
	double W[4],WL[4],WR[4],WS[4];	//�����ʂ̔z��(���x[0]�C���x[1]�C����[2]�C����[3])
	double G[8];			//���̒l�ƃ��ɂ��֐��̒�`
	double S;					//�Ռ��g�̑��x
	double DX;					//bin�̒���
	double XPOS;				//x���W
	double P,U;					//Star Region�ł̈��́C���x�����߂邽�߂̃|�C���^����
	double t;					//����
	int i,t2=0;					//���[�v�p�ϐ�
	
	char filename[50];
	
	//�z��̏�����
	for(i=0;i<4;i++){
		W[i]=0.0;
		WL[i]=0.0;
		WR[i]=0.0;
		WS[i]=0.0;
	}
//	printf("check:������\n");
	
	//Set Up������initial_data.dat����ǂݍ���
/*	double DIAPH;		//�����ڐG�s�A���ʂ̈ʒu
	double DOMLEN;		//Tube�̒���
	double CELLS;		//bin�̌�(�𑜓x)
	double GAMMA;		//��M��
	double TIMEOUT;		//�o�͎���
	double MPA;			//�K�i���萔
	FILE *fp;
	fp= fopen("initial_data.dat","r");
	//�����f�[�^�ǂݍ���
	if((fp == NULL))
	{
		printf("�����f�[�^��ǂݍ��߂܂���ł����D\n");
		return 1;
	}
	fscanf(fp,"%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&DOMLEN,&DIAPH,&CELLS,&GAMMA,&TIMEOUT,&(WL[0]),&(WL[1]),&(WL[2]),&(WR[0]),&(WR[1]),&(WR[2]),&MPA);
		//�ϐ��m�F�p�W���o��
		printf("=================================\n");
		printf("Initial Condition.\n");
		printf("DOMLEN=%lf\nDIAPH=%lf\nCELLS=%d\nGAMMA=%lf\nTIMEOUT=%d\nDL=%lf\nUL=%lf\nPL=%lf\nDR=%lf\nUR=%lf\nPR=%lf\nMPA=%lf\n",DOMLEN,DIAPH,CELLS,GAMMA,TIMEOUT,WL[0],WL[1],WL[2],WR[0],WR[1],WR[2],MPA);
		printf("=================================\n");
	fclose(fp);
*/	
	//�������x�C���x�C���͂̓���
	printf("Left Side 'density' 'velocity' 'pressure'=");
	scanf("%lf%lf%lf",&(WL[0]),&(WL[1]),&(WL[2]));
	//�E�����x�C���x�C���͂̓���
	printf("Right Side 'density' 'velocity' 'pressure'=");
	scanf("%lf%lf%lf",&(WR[0]),&(WR[1]),&(WR[2]));
	
	//�ϐ��m�F�p�W���o��
	printf("===========================================\n");
	printf("Initial Condition.\n");
	printf("(DL,UL,PL)=(%lf,%lf,%lf)\n",WL[0],WL[1],WL[2]);
	printf("(DR,UR,PR)=(%lf,%lf,%lf)\n",WR[0],WR[1],WR[2]);
	printf("DOMLEN=%lf\nDIAPH=%lf\nCELLS=%d\nGAMMA=%lf\nTIMEOUT=%lf\n",DOMLEN,DIAPH,CELLS,GAMMA,TIMEOUT);
	printf("===========================================\n\n\n");
	
	
	//GAMMA�̊֐��̌v�Z
	G[0] = (GAMMA-1.0)/(2.0*GAMMA);
	G[1] = (GAMMA+1.0)/(2.0*GAMMA);
	G[2] = 2.0*GAMMA/(GAMMA-1.0);
	G[3] = 2.0/(GAMMA-1.0);
	G[4] = 2.0/(GAMMA+1.0);
	G[5] = (GAMMA-1.0)/(GAMMA+1.0);
	G[6] = (GAMMA-1.0)/2.0;
	G[7] = GAMMA-1.0; 
	
	//�����̌v�Z
	WL[3] = sqrt(GAMMA*WL[2]/WL[0]);
	WR[3] = sqrt(GAMMA*WR[2]/WR[0]);
	
	//�^��͈���Ȃ��̂Ńv���O�����I��
	if (G[3]*(WL[3]+WR[3]) < (WR[1]-WL[1]))
	{
		//the initial data is such that vacuum is generated.program stopped.
		printf("***Vacuum is generates by data***\n***Program stopped***\n");
		return 1;
	}
	
	//Star Region�ł̈��͂Ƒ��x�����߂�֐��̎��s
	STARPU(&P,&U,WL,WR,WS,G);
	
	//bin���̌v�Z
	DX = DOMLEN/(double)CELLS;
	

	//TIMEOUT���̉������߂�
	while(t<=TIMEOUT){
		FILE *fp;
		sprintf(filename,"exact_%d.dat",t2);
		fp = fopen(filename,"w");
		for(i=1; i<=CELLS; i++){
			XPOS = ((double)i-0.5)*DX;	//x���W�̌v�Z(���ڂ�bin��)
			S = (XPOS-DIAPH)/t;			//����x�ł�TIMEOUT���ɂ����鉹���̌v�Z
			
			//solution at point (x,t)=(XPOS-DIAPH,TIMEPUT) is found
			SAMPLE(S,WS,W,WL,WR,G);
			//exact solution prifiles are written to exact.dat(x���W�Ɩ��x�C���x�C���́C�G�l���M�[�̏o��)
			fprintf(fp, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", XPOS, W[0], W[1], W[2]/MPA, W[2]/(W[0]*G[7]*MPA));
//			printf("FILE OUT:count %d\tXPOS=%lf\n",i,XPOS);
		}
		t = t + dt;
		t2++;
		fclose(fp);
	}
	
	/*�v�����Ԃ̏I��*/
	end=clock();
	printf("endtime:%d[ms]\n",end);
	printf("uptime:%d[ms]\n",end-start);
	
	return 0;
	
}
