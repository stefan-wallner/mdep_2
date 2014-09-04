#include<string>
#include<iostream>


#include "minimize.h"
#include"../chi_squared/breitWigners.h"
#include "yaml-cpp/yaml.h"
#include "currentDateTime.h"

#include"invert33.h"



std::string C="../13w_11t_testload.yaml";


int main(int argc, char* argv[]){
	std::cout<<"--------start------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;
	minimize Chi2(C);

	//int n = Chi2.getParNumber("G_a1(1260)");
	//Chi2.setParLimits(n,.13,13.);
	Chi2.printStatus();
	//return 0;

	int nPar = Chi2.getNpar();
	int nCpl = Chi2.getNanc();
	int nBra = Chi2.getNbra();
	int nTot = Chi2.getNtotAnc();

	std::cout<< "nPar: "<<nPar<<"; nCpl: "<<nCpl<<"; nBra: "<<nBra<<" => nTot: "<<nTot<<std::endl;

	Chi2.conjugate();
//	Chi2.initCouplings();
	print_vector(Chi2.getParameters());


	double params[] = {1.41204, -1.06975, -0.00748951, 0.0380807, 3.3125, -2.89382, -1.54082, 0.803489, 0.00969979, -0.0418916, -3.49037, 1.92446, -0.365443, -1.64996, 0.0421912, 0.00731507, -0.702389, -3.52486, 1.56651, 0.473571, -0.0379001, 0.0522614, 2.9242, 1.18575, -1.37183, 0.801453, -0.0181871, -0.0704837, -2.56476, 0.971844, 1.52324, -0.234316, -0.0175631, 0.0918839, 2.27319, 0.236614, -1.429, -0.3788, 0.0478653, -0.0855039, -1.44503, -1.0808, 0.280363, 1.33411, -0.0942537, -0.0137978, -0.556254, 1.30749, -0.366226, -1.09952, 0.0941395, -0.01343, 0.7383, -0.913639, -0.0014982, -0.679641, 0.0971191, -0.0280913, 1.31206, 0.661233, -0.265434, 0.216036, -0.0225412, -0.0581884, 1.10534, -0.86583, 1.2794, 0.40958, 2.0108, 0.32212, 1.4795, -4.3912, 12.77, -23.516, 1.4093, 0.13925, -5.8146, 2.5961, 1.8005, 0.21337, -3.3995, 1.3136, 0.11089, 1.6704, 0.41407, -0.77708, 0.66388, -1.6812, 0.51046, 2.5439, -0.89233, 1.6521, 0.30451, 1.8335, 0.32157, 2.4491, 2.0204, -3.7016, 3.1433, -2.978, -1.3158, 0.90445, -0.018399, 2.0252, -0.92364, 1.9256, 0.37209, -0.18977, -4.3411, 1.001, 0.00132763, 0.829259, 0.0181491, -15.3911, 2.69214, 119.216, -83.8281, 1460.37, -348.682, -412.197, 6.58593, 167.65, -7.31843, -7.05273, 0.651271, 1546.29, -2450.71, -2127.02, 3600.22, 87.4355, -98.3871, -43.6377, 7.59307, 310.204, 32.4695, 139.365, -40.1092, 9670.56, -4562.28, 909.894, -355.194};

	for (int i=0;i<nTot;i++){
		Chi2.setParameter(i,params[i]);
	};


	std::vector<double> parr(nTot);
	for (int i=0;i<nTot;i++){
		parr[i] = params[i];
	};

	std::cout<<"chi2(): "<<Chi2()<<std::endl;
	std::cout<<"chi2(const double*): "<<Chi2(params)<<std::endl;
	std::cout<<"chi2(std::vector<double>): "<<Chi2(parr)<<std::endl;
	std::cout<<"should be: 46000.3, if everything works (6waves)"<<std::endl;
	std::cout<<"should be: 86023.5 , if everything works (13waves_testload)"<<std::endl;

	Chi2();
	Chi2.write_plots("plots.txt",5);

/*
	std::cout<<"--------init-------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;

	Chi2.initCouplings();
	std::cout<<"--------fit--------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;

	for (int i=2*nCpl;i<nTot;i++){
		std::cout <<"RelPar('"<<Chi2.getParName(i)<<"'):"<<std::endl;
		Chi2.relPar(i);
	};
	std::cout<<Chi2.fit()<<std::endl;
	Chi2.printStatus();
	print_vector(Chi2.getParameters());
	std::cout<<"--------stop-------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;
*/
	return 0;
};

