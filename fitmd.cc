#include<string>
#include<iostream>


#include"minimize.h"
#include"../chi_squared/breitWigners.h"
#include"yaml-cpp/yaml.h"
#include"currentDateTime.h"
#include"matrix_utilities.h"


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
	size_t nPar = Chi2.method()->Waveset()->getNpar();
	size_t nCpl = Chi2.method()->nCpl();
	size_t nBra = Chi2.method()->Waveset()->nBranch();
	size_t nTot = Chi2.method()->nTot();
	size_t nIso = Chi2.method()->Waveset()->getNiso();

	Chi2.printStatus();
	Chi2.method()->Waveset()->printParameters();

	std::cout<< "nPar: "<<nPar<<"; nCpl: "<<nCpl<<"; nBra: "<<nBra<<"; nIso: "<<nIso<<" => nTot: "<<nTot<<std::endl;

//	for(size_t i =0;i<2*nCpl;++i){
//		Chi2.setParameter(i,0.);
//	};
//	print_vector(Chi2.method()->parameters());

//	double mincpl0[] = {-0.0939535, -1.77798, 0.0103546, 0.0239377, -0.551942, -4.35815, -3.18178, -2.0867, 0.992887, 2.85563, 6.02061, 26.0761, 5.73455, -0.62779, -0.981888, -0.455697, 0.595485, -0.176115, 0.143274, 0.510437, -526.927, -922.296, 130.548, 107.55, -73.6623, -122.349, -49.0378, -113.996, 3.49182, 2.94789, -4.71697, -23.504, -9.99731, -13.23, -1.15653, -1.05145, -2.6136, -2.93982, -3212.45, 372.683, 4643.27, -780.66, -202.638, -73.6321, -144.038, -16.5198, 38.3427, 27.3795, 8.35382, -59.6984, -128.853, -322.555, -138.394, -95.5949, -31.9701, 13.0923, -20.6158, -9.41411, 0.70602, -13.6373, -8.73475, -7.51873, -10611.3, -2273.08, -2200.64, 582.229, -912.666, -261.219, -119.908, -1928.88};
//	std::vector<std::complex<double> > ccppll(nCpl,std::complex<double>(0.,0.));
//	for (size_t i=0;i<35;++i){
//		ccppll[i] = std::complex<double>(mincpl0[2*i],mincpl0[2*i+1]);	
//	};
//	print_vector(Chi2.method()->get_branchings(ccppll,std::vector<double>(),std::vector<double>()));

//	for (size_t i=1;i<11;++i){
//		Chi2.method()->Waveset()->setEvalTbin(i,false);
//	};


//	std::cout <<"------->>>>>> "<<(*Chi2.method())()<<std::endl;
//	return 0;	

//	Chi2.findRandRange();
//	return 0;
	std::cout<<"--------init-------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;
	Chi2.initCouplings(2);


	std::cout<<"-------inited------"<<std::endl;
	std::cout<<currentDateTime()<<std::endl;
	std::cout<<"-------------------"<<std::endl;


/*	Chi2.relPar(nTot-1);
	Chi2.relPar(nTot-2);
	Chi2.relPar(nTot-3);
	Chi2.relPar(nTot-4);
	Chi2.method()->write_plots("plots.txt",0);
	Chi2.fit();
	Chi2.method()->write_plots("plots.txt",0);
*/
//	Chi2.open_output("./out_out.out");
//	std::cout<<Chi2()<<std::endl;
//	Chi2.close_output();

//	Chi2.relPar(4);
//	Chi2.relPar(5);
//	Chi2.relPar(6);
//	Chi2.relPar(7);

	return 0;
};

