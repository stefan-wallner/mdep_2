#include<string>
#include<iostream>


#include"minimize.h"
#include"../chi_squared/breitWigners.h"
#include"yaml-cpp/yaml.h"
#include"currentDateTime.h"
#include"AB.h"
#include"invert33.h"
#include"cholesky.h"



//std::string C="../13w_11t_testload.yaml";
//std::string C="/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/chi_squared_retry/6waves_2014-09-05_10:08:53.191171.yaml";
//std::string C="/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/chi_squared_retry/card_test_deiso_2014-09-11_10.38.02.497167.yaml";
int main(int argc, char* argv[]){
	std::string C = std::string(argv[1]);
	std::cout<<"--------start------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;
	minimize Chi2(C);

	//int n = Chi2.getParNumber("G_a1(1260)");
	//Chi2.setParLimits(n,.13,13.);
	//Chi2.printStatus();
	//return 0;
	int nPar = Chi2.method()->Waveset()->getNpar();
	int nCpl = Chi2.method()->getNanc();
	int nBra = Chi2.method()->Waveset()->nBranch();
	int nTot = Chi2.method()->getNtotAnc();
	int nIso = Chi2.method()->Waveset()->getNiso();

	Chi2.printStatus();
	Chi2.method()->Waveset()->printParameters();

	std::cout<< "nPar: "<<nPar<<"; nCpl: "<<nCpl<<"; nBra: "<<nBra<<"; nIso: "<<nIso<<" => nTot: "<<nTot<<std::endl;
	
	std::cout<<"--------init-------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;
	Chi2.initCouplings();


	std::cout<<"-------inited------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;
	Chi2.relPar(nTot-1);
	Chi2.relPar(nTot-2);
	Chi2.relPar(nTot-3);
	Chi2.relPar(nTot-4);
	Chi2.method()->write_plots("plots.txt",0);
	Chi2.fit();
	Chi2.method()->write_plots("plots.txt",0);
//	Chi2.open_output("./out_out.out");
//	std::cout<<Chi2()<<std::endl;
//	Chi2.close_output();

//	Chi2.relPar(4);
//	Chi2.relPar(5);
//	Chi2.relPar(6);
//	Chi2.relPar(7);

	return 0;
};

