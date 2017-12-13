#include "riemann_function.hpp"

void GUESSP(quantity *WL, quantity *WR, quantity *WS, double G[8]){
	//変数定義
	double PPV;					//pressure based on primitive variables
	double PMIN,PMAX;		//左右の圧力を比較し，大小を比較するためのもの
	double QMAX;				//圧力比
	double QUSER = 2.0;	//スイッチングパラメータ
	double PQ,PTL,PTR;	//two rarefaction approximationの際に使うもの
	double GEL,GER;			//two shock appriximationの際に使うもの

	PPV = ((WL->get_pressure()+WR->get_pressure())*0.5)
			+ ((WL->get_velocity()-WR->get_velocity())
				*(WL->get_density()+WR->get_density())
				*(WL->get_sound()+WR->get_sound())*0.125);
	PPV = fmax(0.0, PPV);
	PMIN = fmin(WL->get_pressure(), WR->get_pressure());
	PMAX = fmax(WL->get_pressure(), WR->get_pressure());
	QMAX = PMAX/PMIN;

	if(QMAX<QUSER && (PMIN<PPV && PPV<PMAX)){
		WS->set_pressure(PPV);
	}else if(PPV < PMIN){
		//Two Rarefanction approximation
		PQ = pow((WL->get_pressure()/WR->get_pressure()),G[0]);
		WS->set_velocity((PQ*WL->get_mach() + WR->get_mach() + G[3]*(PQ-1.0))
											/ (PQ/WL->get_sound()+1.0/WR->get_sound()));
		PTL = 1.0+G[6]*(WL->get_velocity()-WS->get_velocity())/WL->get_sound();
		PTR = 1.0+G[6]*(WS->get_velocity()-WR->get_velocity())/WR->get_sound();
		WS->set_pressure(0.5*(pow(WL->get_pressure()*PTL,G[2])+pow(WR->get_pressure()*PTR,G[2])));
		cout << "Two-Rarefaction approximation selected." << endl;
	}else{
		//select two-shock Riemann solver with PVRS
		GEL = sqrt((G[4]/WL->get_density())/(G[5]*WL->get_pressure()+PPV));
		GER = sqrt((G[4]/WR->get_density())/(G[5]*WR->get_pressure()+PPV));
		WS->set_pressure((GEL*WL->get_pressure()+GER*WR->get_pressure()-(WR->get_velocity()-WL->get_velocity()))/(GEL+GER));
		cout << "Two-Shock approximation selected." << endl;
	}
}

void PREFUNL(double *FL, double *FLD, quantity *WS, quantity *WL, double G[8]){
	double AL,BL,PRATL,QRTL;
	if(WS->get_pressure() < WL->get_pressure()){
		//rarefaction wave
		PRATL = WS->get_pressure()/WL->get_pressure();
		*FL = G[3]*WL->get_sound()*(pow(PRATL,G[0]) - 1.0);
		*FLD = (1.0/(WL->get_density()*WL->get_sound()))/pow(PRATL,G[1]);
	}else{
		//shock wave
		AL = G[4]/WL->get_density();
		BL = G[5]*WL->get_pressure();
		QRTL = sqrt(AL/(BL + WS->get_pressure()));
		*FL = (WS->get_pressure() - WL->get_pressure())*QRTL;
		*FLD = (1.0 - (WS->get_pressure() - WL->get_pressure())/(2.0*(BL + WS->get_pressure())))*QRTL;
	}
}

void PREFUNR(double *FR, double *FRD, quantity *WS, quantity *WR, double G[8]){
	double AR,BR,PRATR,QRTR;
	if (WS->get_pressure() < WR->get_pressure()){
		//rarefaction wave
		PRATR = WS->get_pressure()/WR->get_pressure();
		*FR= G[3]*WR->get_sound()*(pow(PRATR,G[0]) - 1.0);
		*FRD = (1.0/(WR->get_density()*WR->get_sound()))/pow(PRATR,G[1]);
	}else{
		//shock wave
		AR = G[4]/WR->get_density();
		BR = G[5]*WR->get_pressure();
		QRTR = sqrt(AR/(BR + WS->get_pressure()));
		*FR = (WS->get_pressure() - WR->get_pressure())*QRTR;
		*FRD = (1.0 - (WS->get_pressure() - WR->get_pressure())/(2.0*(BR + WS->get_pressure())))*QRTR;
	}
}

void STARPU(double *P, double *U, quantity *WL, quantity *WR, quantity *WS, double G[8]){
	int i, NRITER = 10;								//反復回数
	double TOLPRE = 1/pow(10.0,6.0);	//p0であるという条件
	double UDIFF;											//左右波の速度差
	double FL,FLD,FR,FRD;							//左右の圧力に関する補助関数
	double CHANGE;										//圧力の差分

	//Guessed value PSTART is computed
	GUESSP(WL, WR ,WS, G);
	cout << "guess result PM = " << WS->get_pressure() << endl;
	UDIFF = WR->get_velocity() - WL->get_velocity();

	cout << "===========================================" << endl;
	cout << "Interation number:Change" << endl;
	cout << "===========================================" << endl;
	for(i=1; i<=NRITER; i++){
		PREFUNL(&FL, &FLD, WS, WL, G);
		PREFUNR(&FR, &FRD, WS, WR, G);
		*P = WS->get_pressure()-(FL+FR+UDIFF)/(FLD+FRD);
		CHANGE = 2.0*fabs((*P-WS->get_pressure()))/(*P+WS->get_pressure());
		cout << "i = " << i << "CHANGE = " << CHANGE << endl;
		if(CHANGE<TOLPRE)	continue;
		if(*P<0.0) *P=TOLPRE;
		WS->set_pressure((*P));
	}
	printf("Divergence in Newton-Raphson interation\n");
	*U = (WL->get_velocity()+WR->get_velocity()+FR-FL)*0.5;
	WS->set_velocity((*U));
	//Compute velocity in star region
	cout << "\n===========================================" << endl;
	cout << "Pressure\t Velocity" << endl;
	cout << "===========================================" << endl;
	cout << *P/MPA << " " << *U << endl;
	cout << "===========================================" << endl;
}

void SAMPLE(double S, quantity *WS, quantity *W, quantity *WL, quantity *WR, double G[4]){
	double SHL,STL;		//Left rarefaction head,tail
	double SHR,STR;		//Right rarefaction head,tail
	double CML,CMR;		//rarefaction内の音速の計算
	double PML,PMR;		//SL,SRを計算するためのもの
	double SL,SR;			//左右のshockの速度

	if (S < WS->get_velocity()){
		//Sampling point lies to the left of the contact discontinuity
		if(WS->get_pressure() < WL->get_pressure()){
			//left rarefaction
			SHL=WL->get_velocity()-WL->get_sound();	//左側rarefactionの速度
			if(S < SHL){
				//samples point is lecft data state
				W->set_density(WL->get_density());
				W->set_velocity(WL->get_velocity());
				W->set_pressure(WL->get_pressure());
			}else{
				CML = WL->get_sound()*pow((WS->get_pressure()/WL->get_pressure()),G[0]);
				STL = WS->get_velocity() - CML;
				if(S > STL){
					//samples point is star left state
					W->set_density(WL->get_density()*pow((WS->get_pressure()/WL->get_pressure()),(1.0/GAMMA)));
					W->set_velocity(WS->get_velocity());
					W->set_pressure(WS->get_pressure());
				}else{
					//sampled point is inside left fan
					W->set_velocity(G[4]*(WL->get_sound() + G[6]*WL->get_velocity() + S));
					W->set_sound(G[4]*(WL->get_sound() + G[6]*(WL->get_velocity() - S)));
					W->set_density(WL->get_density()*pow((W->get_sound()/WL->get_sound()),G[3]));
					W->set_pressure(WL->get_pressure()*pow((W->get_sound()/WL->get_sound()),G[2]));
				}
			}
		}else{
			//left shock
			PML = WS->get_pressure()/WL->get_pressure();
			SL = WL->get_velocity()-WL->get_sound()*sqrt(G[1]*PML+G[0]);
			if(S < SL){
				//sampled point is left data state
				W->set_density(WL->get_density());
				W->set_velocity(WL->get_velocity());
				W->set_pressure(WL->get_pressure());
			}else{
				//sampled point is star left state
				W->set_density(WL->get_density()*(PML + G[5])/(PML*G[5] + 1.0));
				W->set_velocity(WS->get_velocity());
				W->set_pressure(WS->get_pressure());
			}
		}
	}else if(WS->get_pressure() > WR->get_pressure()){
		PMR = WS->get_pressure()/WR->get_pressure();
		SR = WR->get_velocity() + WR->get_sound()*sqrt(G[1]*PMR + G[0]);
		if(S > SR){
			//Sampled point is right data state
			W->set_density(WR->get_density());
			W->set_velocity(WR->get_velocity());
			W->set_pressure(WR->get_pressure());
		}else{
			//Sampled point is Star Right state
			W->set_density(WR->get_density()*(PMR + G[5])/(PMR*G[5] + 1.0));
			W->set_velocity(WS->get_velocity());
			W->set_pressure(WS->get_pressure());
		}
	}else{
		//right rarefaction
		SHR = WR->get_velocity()+WR->get_sound();	//head wave
		if(S > SHR){
			//sampled point is right data state
			W->set_density(WR->get_density());
			W->set_velocity(WR->get_velocity());
			W->set_pressure(WR->get_pressure());
		}else{
			CMR = WR->get_sound()*pow((WS->get_pressure()/WR->get_pressure()),G[0]);
			STR = WS->get_velocity() + CMR;
			if(S < STR){
				//sampled point is star right state
				W->set_density(WR->get_density()*pow((WS->get_pressure()/WR->get_pressure()),(1.0/GAMMA)));
				W->set_velocity(WS->get_velocity());
				W->set_pressure(WS->get_pressure());
			}else{
				//Sampled point is inside right fan
				W->set_velocity(G[4]*(-WR->get_sound() + G[6]*WR->get_velocity() + S));
				W->set_sound(G[4]*(WR->get_sound() - G[6]*(WR->get_velocity() - S)));
				W->set_density(WR->get_density())*pow((W->get_sound()/WR->get_sound()),G[3]);
				W->set_pressure(WR->get_pressure()*pow((W->get_sound()/WR->get_sound()),G[2]));
			}
		}
	}
}
