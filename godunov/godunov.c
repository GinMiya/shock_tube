#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//#define GAMMA	1.4		//��M��
#define GAMMA	(5.0/3.0)	//��M��
#define DOMLEN	1.0		//tube�̒���
#define DIAPH	0.5		//�d�؂�̈ʒu
#define CELLS	100		//�Z���̐�
#define TEND	0.002	//�I������
#define MPA		1.0		//�K�i���萔
#define CFLCOE	0.9		//CFL�W��

//�֐���`
//�萔�o��
void PRINT(){
	//�l�̕W���o��
	printf("=======================================\n");
	printf("Dates of This Simulation\n");
	printf("=======================================\n");
	printf("CFLCOE=%lf\n",CFLCOE);
	printf("DOMLEN=%lf\n",DOMLEN);
	printf("DIAPH=%lf\n",DIAPH);
	printf("CELLS=%d\n",CELLS);
	printf("GAMMA=%lf\n",GAMMA);
}
//GAMMA�Ɋւ���֐�
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
//�^��̏�������
int VACUUM(double INIWL[4], double INIWR[4], double G[8]){
	//�^��͈���Ȃ��̂Ńv���O�����I��
	if(G[3]*(INIWL[3]+INIWR[3])<(INIWR[1]-INIWL[1])){
		//the initial data is such that vacuum is generated.program stopped.
		printf("***Vacuum is generates by data***\n***Program stopped***\n");
		return 1;
	}//end if
}//end VACUUM
//���������ʂ̏o��
void INICOND(double INIWL[4], double INIWR[4], double G[8]){
	//�����̌v�Z�Əo�́�CFL�����Ɏg��
	INIWL[3] = sqrt(GAMMA*INIWL[2]/INIWL[0]);
	INIWR[3] = sqrt(GAMMA*INIWR[2]/INIWR[0]);
	printf("Initial Condition\n");
	printf("(DL,UL,PL,AL)=(%lf,%lf,%lf,%lf)\n",INIWL[0],INIWL[1],INIWL[2],INIWL[3]);
	printf("(DR,UR,PR,AR)=(%lf,%lf,%lf,%lf)\n",INIWR[0],INIWR[1],INIWR[2],INIWR[3]);
	printf("=============================================================================\n\n");
}
//������Ԃ̐ݒ�֐�CELL�ɕ����ʂ�����̂Ə�����/
void INITIAL(double U[][3], double INIWL[4], double INIWR[4], double DX, double G[8]){
	int i,j,k;
	//�ۑ��ʂ̏�����
	for(i=0;i<CELLS+2;i++){
		for(k=0;k<3;k++){
			U[i][k]=0.0;
		}//end for k
	}//end for i
	//���E�����̕����ʂ�����(�ۑ��`)/�����͔����̃Z�������E�ɂ��ĕ����ʂ̕��z(�Ռ��g�ǖ��̂���)
	for(i=1;i<DIAPH/DX+1;i++){
		//���x
		U[i][0]=INIWL[0];
		//���x�~���x
		U[i][1]=INIWL[0]*INIWL[1];
		//�P�ʎ��ʂ�����̑S�G�l���M�[
		U[i][2]=INIWL[2]/G[7]+0.5*INIWL[0]*INIWL[1]*INIWL[1];
	}//end for i
	for(i=DIAPH/DX+1;i<CELLS+1;i++){
		//���x
		U[i][0]=INIWR[0];
		//���x�~���x
		U[i][1]=INIWR[0]*INIWR[1];
		//�P�ʎ��ʂ�����̑S�G�l���M�[
		U[i][2]=INIWR[2]/G[7]+0.5*INIWR[0]*INIWR[1]*INIWR[1];
	}//end for i
}//end INITIAL
//���E�����֐�
void BOUNDI(double U[][3]){
	int k;
	for(k=0;k<3;k++){
		U[0][k]=0.0;
		U[CELLS+1][k]=0.0;
		U[0][k]=U[1][k];
		U[CELLS+1][k]=U[CELLS][k];
	}//end for
}//end BOUNDI
//CFL����
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
*�������牺��Riemann���ł��� *
********************************/
void GUESSP(double WL[4], double WR[4], double WS[4], double G[8]){
	//�ϐ���`
	double PPV;			//pressure based on primitive variables(4.47)���Q��
	double PMIN,PMAX;	//���E�̈��͂��r���C�召���r���邽�߂̂���
	double QMAX;		//���͔�
	double QUSER = 2.0;	//�X�C�b�`���O�p�����[�^
	double PQ,PTL,PTR;	//two rarefaction approximation�̍ۂɎg������
	double GEL,GER;		//two shock appriximation�̍ۂɎg������
	//compute guess pressure from PVRS Riemann solver
	PPV = (0.5*(WL[2]+WR[2])/2.0)+(0.125*(WL[1]-WR[1])*(WL[0]+WR[0])*(WL[3]+WR[3]));	//PPV�̌v�Z(4.47)��
	PPV = fmax(0.0, PPV);		//0.0��PPV�̑傫���ق����Ƃ�
	PMIN = fmin(WL[2], WR[2]);	//���E�̏������ق����Ƃ�
	PMAX = fmax(WL[2], WR[2]);	//���E�̑傫���ق����Ƃ�
	QMAX = PMAX/PMIN;			//���͔�
	//select two-shock Riemann solver with PVRS as estimate (4.48)�� ���ꂪ�̗p�����͂�
		GEL = sqrt((G[4]/WL[0])/(G[5]*WL[2]+PPV));
		GER = sqrt((G[4]/WR[0])/(G[5]*WR[2]+PPV));
		WS[2] = (GEL*WL[2]+GER*WR[2]-(WR[1]-WL[1]))/(GEL+GER);
}//end GUESSP
//Riemann Solver�ł�f�֐��̌v�Z(����)
void PREFUNL(double *FL, double *FLD, double WS[4], double WL[4], double G[8]){
	double AL,BL,PRATL,QRTL;	//�萔A_K,B_K��F�̂��߂̂���
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
//Riemann Solver�ł�f�֐��̌v�Z(�E��)
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
//Star Region�ł̈��͂Ƒ��x�����߂邽�߂̊֐�
void STARPU(double WL[4], double WR[4], double WS[4], double G[8]){
	int j;
	int NRITER = 10;					//������
	double TOLPRE = 1.0/pow(10.0,6.0);	//p0�ł���Ƃ�������
	double UDIFF;						//���E�g�̑��x��
	double FL,FLD,FR,FRD;				//���E�̈��͂Ɋւ���֐�
	double CHANGE;						//���͂̍���
	double PS,US;						//StarRegion�̈��͂Ƒ��x�̐���l
	//Guessed value PSTART is computed
	GUESSP(WL,WR,WS,G);				//�����񐔂�����悤�Ȉ��͂̐��@
	UDIFF = WR[1] - WL[1];			//�g�̍��E�ł̑��x���v�Z
	for(j=1; j<=NRITER; j++){
		PREFUNL(&FL, &FLD, WS, WL, G);				//�֐�F�̌v�Z(����)
		PREFUNR(&FR, &FRD, WS, WR, G);				//�֐�F�̌v�Z(�E��)
		PS = WS[2]-(FL+FR+UDIFF)/(FLD+FRD);			//�V���Ȉ��͂̌v�Z(Newton�@)
		CHANGE = 2.0*fabs((PS-WS[2]))/(PS+WS[2]);	//���͂̕ω��ʂ��v�Z
		if(CHANGE<TOLPRE){
			continue;		//�ω��ʂ�10^-6�����������Ƃ���interation���Ȃ�
		}
		if(PS<0.0){
			PS=TOLPRE;		//���͂����ɂȂ����炻�̂Ƃ��͍ł�������10^-6�Ƃ���
		}
		WS[2]=PS;
	}						//interation�̏I��
	//interation�ł킩����star region�ł̈��͂��g���đ��x���v�Z
	US = (WL[1] + WR[1] + FR - FL)/2.0;
	WS[1]=US;
}//end STARPU
//�Ǐ�Riemann���
void RIEMANN(double WS[4], double W[4], double WL[4], double WR[4], double G[8], double XPOS){
	double SHL,STL;		//Left rarefaction head,tail
	double SHR,STR;		//Right rarefaction head,tail
	double CML,CMR;		//rarefaction���̉����̌v�Z
	double PML,PMR;		//SL,SR���v�Z���邽�߂̂���
	double SL,SR;		//���E��shock�̑��x
	//printf("RIEMANN pressure velocity in StarRegion=%lf\t%lf\n",WS[2],WS[1]);
	if(WS[1]>0.0){
		if(WS[2]>WL[2]){
			//shock wave
			PML=WS[2]/WL[2];
			//����shock���x
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
			SHL=WL[1]-WL[3];	//����rarefaction�̑��x
			CML=WL[3]*pow((WS[2]/WL[2]),G[0]);//a*L�̌v�Z(4.55)
			STL=WS[1]-CML;	//����rarefaction�̑��x
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
			SHR=WR[1]+WR[3];	//�E��rarefaction head�̑��x
			CMR=WR[3]*pow((WS[2]/WR[2]),G[0]);	//a*R�̌v�Z(4.62)
			STR=WS[1]+CMR;		//�E��rarefaction tail�̑��x
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
//�t���b�N�X�v�Z
void FLUXES(double U[][3], double FLUX[][3], double WL[4], double WR[4], double W[4], double WS[4], double G[8], double DX, double XPOS){
	int i,j,k;
	for(i=0;i<CELLS+1;i++){
		//Step4-1:FLUX�̏�����
		for(k=0;k<3;k++){
			W[k]=0.0;
			WS[k]=0.0;
			FLUX[i][k]=0.0;
		}//end for
		//Step4-2:���W���Ƃ�Rieamann����������Flux�����߂�
		XPOS=(i+0.5)*DX;
		//�ۑ��ʂ𑬓x�C���͂ɂ��ď����Ȃ����Ă��K�v������
			//Riemann���p�������x,���x,����,����
			WL[0]=U[i][0];
			WL[1]=U[i][1]/(U[i][0]);
			WL[2]=(U[i][2]-0.5*WL[0]*WL[1]*WL[1])*G[7];
			WL[3]=sqrt(fabs(GAMMA*WL[2]/WL[0]));
			//Riemann���p�E�����x,���x,����,����
			WR[0]=U[i+1][0];
			WR[1]=U[i+1][1]/(U[i+1][0]);
			WR[2]=(U[i+1][2]-0.5*WR[0]*WR[1]*WR[1])*G[7];
			WR[3]=sqrt(fabs(GAMMA*WR[2]/WR[0]));
		//Riemann���������ɂ������Ă�StarRegion�ł̈���,���x�̌��ς���
		STARPU(WL,WR,WS,G);
		//RIEMANN�̎��s�Ŗ��x,���x,���͂����Ƃ߂�
		RIEMANN(WS,W,WL,WR,G,XPOS);
		//Eular�������̃t���b�N�X
		FLUX[i][0]=W[0]*W[1];
		FLUX[i][1]=W[0]*W[1]*W[1]+W[2];
		FLUX[i][2]=W[1]*0.5*(G[2]*W[2]+W[0]*W[1]*W[1]);
	}//end for
}//end FLUX
//���̃^�C���X�e�b�v�փA�b�v�f�[�g
void UPDATE(double DTODX, double U[][3], double FLUX[][3]){
	int i,j,k;
	for(i=1;i<CELLS+1;i++){
		U[i][0] = U[i][0]+DTODX*(FLUX[i-1][0]-FLUX[i][0]);
		U[i][1] = U[i][1]+DTODX*(FLUX[i-1][1]-FLUX[i][1]);
		U[i][2] = U[i][2]+DTODX*(FLUX[i-1][2]-FLUX[i][2]);
	}//end for
}//end UPDATE

//main�֐�
int main(void){
	//�萔��`
	int i,j,k,N=1;
	double G[8];						//GAMMA�֐�
	double INIWL[4],INIWR[4];			//���������ʔz��
	double U[CELLS+2][3],FLUX[CELLS+2][3];	//�ۑ��ʂƃt���b�N�X
	double W[4],WS[4];					//Riemann���̉���StarRegion�̕�����
	double WL[4],WR[4];					//Riemann���̍��E�̕�����
	double FL,FLD,FR,FRD;				//Riemann���p�̈��͊֐�
	double r,u,p,e;						//�o�͗p���x,���x,����,�����G�l���M�[
	double P_,U_;						//STARPU�̂��߂̃|�C���^
	double DX=0.0,DT=0.0,DTODX=0.0;		//�Z����,�^�C���X�e�b�v,DT/DX
	double XPOS;						//x���̈ʒu
	double TIME=0.0;					//����
	//�t�@�C���o�͏���
	char dataname[50];
	FILE *fp;FILE *fp0;
	fp0=fopen("godunov_0.dat","w");
	//GAMMA�Ɋւ���֐��̒�`
	GFUNC(G);
	//�e�萔�̏o��
	PRINT();
	//���������ʂ̓���
		printf("=============================================================================\n");
		//�������x�C���x�C���͂̓���
		printf("Left  Side 'density' 'velocity' 'pressure'=");
		scanf("%lf%lf%lf",&(INIWL[0]),&(INIWL[1]),&(INIWL[2]));
		//�E�����x�C���x�C���͂̓���
		printf("Right Side 'density' 'velocity' 'pressure'=");
		scanf("%lf%lf%lf",&(INIWR[0]),&(INIWR[1]),&(INIWR[2]));
		printf("=============================================================================\n");
	//�^�󔻒�
	VACUUM(INIWL,INIWR,G);
	//���������ʂ̊m�F
	INICOND(INIWL,INIWR,G);
	//�Z���̒���/�����ۑ��ʂ̊i�[
	DX=DOMLEN/((double)CELLS);
	INITIAL(U,INIWL,INIWR,DX,G);
	//���������ʂ̏o��
	for(i=1;i<CELLS+1;i++){
		//�ۑ��ʂ𕨗��ʂŏ���
		r=U[i][0];
		u=U[i][1]/r;
		p=(U[i][2]-0.5*r*u*u)*G[7];
		e=p/(r*G[7]*MPA);
		fprintf(fp0,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",FLUX[i][0],FLUX[i][1],FLUX[i][2],U[i][0],U[i][1],U[i][2]);
		//fprintf(fp0,"%lf\t%lf\t%lf\t%lf\t%lf\n",i*DX,r,u,p,e);
	}//end for
	fclose(fp0);
	//���ԃ��[�v
	do{
		//���E����
		BOUNDI(U);
		////CFL��������^�C���X�e�b�v����
		CFLCOND(U,&DT,DX,G);
		DTODX=DT/DX;
		////���Ԃ�DT�����i�߂�
		TIME+=DT;
		////�����ʂ��v�Z,�o��
		printf("=============================================================================\n");
		printf("TIME=%lf\tNo.%d TIMESTEP DT=%lf\t DT/DX=%lf\n",TIME,N,DT,DTODX);
		printf("-----------------------------------------------------------------------------\n");
		////�t���b�N�X�̌v�Z
		FLUXES(U,FLUX,WL,WR,W,WS,G,DX,XPOS);
		////�ۑ��ʂ̃A�b�v�f�[�g
		UPDATE(DTODX,U,FLUX);
		sprintf(dataname,"godunov_%lf.dat",TIME);
		fp=fopen(dataname,"w");
		for(i=1;i<CELLS+1;i++){
			//�ۑ��ʂ𕨗��ʂŏ���
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
