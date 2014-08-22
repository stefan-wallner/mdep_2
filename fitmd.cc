#include<string>
#include<iostream>
#include<map>


#include "minimize.h"
#include"../chi_squared/breitWigners.h"
#include "yaml-cpp/yaml.h"

#include"invert33.h"

double MIN_STEP_SIZE = 0.00001;

const char* wave_definitions = "../wave_specifications.yaml";
const char* parametrizations = "../parametrizations_13_waves.yaml";
//const char* parametrizations = "../parametrizations.yaml";//Commented version with 6 waves works

minimize init_chi2(std::string card){
	YAML::Node param  = YAML::LoadFile(parametrizations);	
	YAML::Node defs   = YAML::LoadFile(wave_definitions);
	YAML::Node waveset= YAML::LoadFile(card);
	minimize Chi2 = minimize();

// Set global phase space paramerizations
	bool sgps = false;
	if (waveset["standard_phase_space"]){
		if (waveset["standard_phase_space"].as<bool>()){
			Chi2.setGlobalPhaseSpace(10);
			sgps=true;
		};
	};
	double min_step_size;
	if (waveset["min_step_size"]){
		min_step_size = waveset["min_step_size"].as<double>();
	}else{
		min_step_size = MIN_STEP_SIZE;
	};

	int nPS = waveset["global_phase_space"].size();
	for (int i=0;i<nPS;i++){
		Chi2.setGlobalPhaseSpace(waveset["global_phase_space"][i].as<int>());
		sgps=true;
	};
	if (not sgps){
		std::cout<<"Warning: No global phase space set" << std::endl;
	};
// Define functions
	std::map<std::string,int> fMap;
	int fCount = 1; // Starts as one, since map[nonExistingKey] = 0 -> use this to determine, whether a function exists.
	int pCount = 0;
	int cCount = 0;
	int nWaves = waveset["waves"].size();
	std::vector<double> pVals;
	std::vector<double> pSteps;
	std::vector<bool> pRels;
	for (int i=0;i<nWaves;i++){
		int nFunc = waveset["waves"][i]["parametrizations"].size();
		for (int j=0;j<nFunc;j++){
			std::string fName = waveset["waves"][i]["parametrizations"][j].as<std::string>();
			if(param[fName]){
				if (fMap[fName]==0){
					fMap[fName]=fCount;
					int nParam = param[fName]["parametrization"].as<int>();
					int nPar = getNparsNonConst(nParam);
					int nCon = getNpars(nParam)-getNparsNonConst(nParam);
					bool isTdep = false;
					if (param[fName]["t_dependent"]){
						if (param[fName]["t_dependent"].as<bool>()){
							isTdep = true;
						};
					};
					Chi2.add_func(nParam,isTdep);
					Chi2.setFunctionName(fCount-1,fName);
					fCount++;
					int nParDef = param[fName]["parameters"].size();
					if (nParDef == nPar){
						for(int par=0;par<nPar;par++){
							std::string pName = param[fName]["parameters"][par]["name"].as<std::string>();
							double value = param[fName]["parameters"][par]["value"].as<double>();
							bool rel = false;
							if(param[fName]["parameters"][par]["released"]){
								if(param[fName]["parameters"][par]["released"].as<bool>()){
									rel = true;
								};
							};
							Chi2.setParameterName(pCount,pName);
							if(param[fName]["parameters"][par]["upper_limit"] and param[fName]["parameters"][par]["lower_limit"]){
								double upper = param[fName]["parameters"][par]["upper_limit"].as<double>();
								double lower = param[fName]["parameters"][par]["lower_limit"].as<double>();
								Chi2.setParLimits(par,upper,lower);
							};
							pCount++;
							pVals.push_back(value);
							double step = fabs(0.001*value);
							if (step < min_step_size){
								step = min_step_size;
							};
							pRels.push_back(rel);
							pSteps.push_back(step);
						};
					}else{
						std::cerr<<"Error: Number of defined parameters does not match required number for "<<fName<<std::endl;
					};
					int nConDef = param[fName]["constants"].size();
					if(nConDef == nCon or (nConDef == nCon-1 and isTdep)){
						for (int con=0;con<nCon;con++){
							std::string cName = param[fName]["constants"][con]["name"].as<std::string>();
							Chi2.setConstantName(cCount ,cName);
							if (con < nCon -1 or not isTdep){
								double value  = param[fName]["constants"][con]["value"].as<double>();
								Chi2.setConst(cCount,value);
							};
							cCount++;
						};
					}else{
						std::cerr<<"Error: Number of defined constants does not match required number for "<<fName<<std::endl;
						std::cerr<<nConDef<<" "<<nCon<<std::endl;
					};					
				};
			}else{
				std::cerr << "Error: '"<<fName<<"' not defined in "<<parametrizations<<std::endl;
			};
		};

	};
// Define waves
	double mmin = waveset["waves"][0]["mmin"].as<double>();
	double mmax = waveset["waves"][0]["mmax"].as<double>();
	for (int i=0;i<nWaves;i++){
		Chi2.add_wave();
		std::string wName = waveset["waves"][i]["name"].as<std::string>();
		Chi2.setWaveName(i,wName);
		if (not defs[wName]){
			std::cerr<<"Error: "<< wName<<" not defined in "<<wave_definitions<<std::endl;
		};
		Chi2.setWaveSpin(i,defs[wName]["spin"].as<int>());
		double mmin_act = waveset["waves"][i]["mmin"].as<double>();
		double mmax_act = waveset["waves"][i]["mmax"].as<double>();
		if (mmin_act<mmin){
			std::cerr << "Error: Lower limit of '"<<waveset["waves"][i]["name"]<<"' below anchor mass limit"<<std::endl;
		};
		if (mmax_act>mmax){
			std::cerr << "Error: Upper limit of '"<<waveset["waves"][i]["name"]<<"' above anchor mass limit"<<std::endl;
		};
		Chi2.setWaveLimits(i,mmin_act,mmax_act);
		Chi2.setWavePhaseSpace(i,defs[wName]["phase_space"].as<int>());
		int nFunc = waveset["waves"][i]["parametrizations"].size();
		for (int func=0;func<nFunc;func++){
			std::string fName = waveset["waves"][i]["parametrizations"][func].as<std::string>();
			int nf = fMap[fName]-1;
			Chi2.add_func_to_wave(i,nf);
		};
	};
// Set up branchings
	std::vector<int> fixed_branchings;
	int nBr=0;
	int nBranch = waveset["branchings"].size();
	for (int i=0;i<nBranch;i++){
		std::vector<int> wave_numbers;
		int nnn = waveset["branchings"][i].size();
		for (int j=0;j<nnn;j++){
			int func_count = 0;
			std::string wnam_b = waveset["branchings"][i][j][0].as<std::string>();
			std::string fnam_b = waveset["branchings"][i][j][1].as<std::string>();
			for (int wave=0;wave<nWaves;wave++){
				std::string wnam_d = waveset["waves"][wave]["name"].as<std::string>();
				int nfunc = waveset["waves"][wave]["parametrizations"].size();
				for (int func=0;func<nfunc;func++){
					std::string fnam_d = waveset["waves"][wave]["parametrizations"][func].as<std::string>();
					if (fnam_d == fnam_b and wnam_d == wnam_b){
						wave_numbers.push_back(func_count);
						};
					func_count++;
				};
			};
		};
		if (wave_numbers.size()<2){
			std::cerr<<"Error: Branching definition with less than two waves encountered."<<std::endl;
		}else{
			for (int cplto = 1; cplto < wave_numbers.size(); cplto++){
				Chi2.couple_funcs(wave_numbers[0],wave_numbers[cplto]);
			};
		};
	};
//Set up data
	double M = waveset["mmin"].as<double>();
	double m_max = waveset["mmax"].as<double>();
	double binwidth = waveset["binwidth"].as<double>();
	std::vector<double> binning;
	binning.push_back(M);
	while(M<m_max){
		M+=binwidth;
		binning.push_back(M);
	};
	print_vector(binning);
	Chi2.setBinning(binning);
	std::vector<double> tbinning;
	int nTbin = waveset["t_binning"].size();
	for (int i=0; i<nTbin;i++){
		tbinning.push_back(waveset["t_binning"][i].as<double>());
	};
	Chi2.setTbinning(tbinning);
	if (waveset["data_files"].size() != nTbin -1){
		std::cerr<<"Error: Number of data-files does not match number of t' bins"<<std::endl;
	};
	if (waveset["coma_files"].size() != nTbin-1){
		std::cerr<<"Error: Number of coma-files does not match number of t' bins"<<std::endl;
	};

	Chi2.update_definitions(); // Has to be called once with all functions set, to get internal handling of parameter numbers right

	for(int i=0;i<nTbin-1;i++){
		Chi2.loadData(i,waveset["data_files"][i].as<std::string>().c_str());
		Chi2.loadComa(i,waveset["coma_files"][i].as<std::string>().c_str());
	};
	int nPar = Chi2.getNpar();
	int nCpl = Chi2.getNanc();
	if (nPar != pVals.size()){
		std::cerr<<"Error: Number of parameters does not match"<<std::endl;
	};
	if (waveset["real_anchor_cpl"]){
		if(waveset["real_anchor_cpl"].as<bool>()){
			Chi2.setParameter(1,0.);
			Chi2.fixPar(1);
		};
	};
	for (int par = 0;par<nPar;par++){
		Chi2.setParameter(2*nCpl+par,pVals[par]);
		Chi2.setStepSize(2*nCpl+par,pSteps[par]);
	};
	std::vector<int> first_branch = Chi2.getFirstBranch();  //  Fix the first branching for each coupling, to remove ambiguities
								// Otherwise B_i C_t would be all the same after B_i -> B_i*D and C_t -> C_t/D
								// For any complex D
	for (int i=0; i<first_branch.size();i++){
		Chi2.setParameter(2*nCpl+nPar+2*first_branch[i]  ,1.); // Re(Br)
		Chi2.setParameter(2*nCpl+nPar+2*first_branch[i]+1,0.); // Im(Br)
		Chi2.fixPar(2*nCpl+nPar+2*first_branch[i]  );
		Chi2.fixPar(2*nCpl+nPar+2*first_branch[i]+1);
	};
	Chi2.setRandomCpl();
	Chi2.setRandomBra();
	Chi2.initialize();
	return Chi2;
};

int main(int argc, char* argv[]){


//	minimize Chi2 = init_chi2("../6w_11t.yaml"); //Commented version with 6 waves works
	minimize Chi2 = init_chi2("../13w_11t.yaml");
	int nPar = Chi2.getNpar();
	int nCpl = Chi2.getNanc();
	int nBra = Chi2.getNbra();
	int nTot = Chi2.getNtotAnc();

	Chi2.printParameters();
//////////////////////////////////////////////////////////
	Chi2.setParameter("M_a1(1260)",1.2794);
	Chi2.setParameter("G_a1(1260)",0.40958);
	Chi2.setParameter("M_a1(1930)",2.0108);
	Chi2.setParameter("G_a1(1930)",0.32212);
	Chi2.setParameter("M_a1(1420)",1.4093);
	Chi2.setParameter("G_a1(1420)",0.13925);
	Chi2.setParameter("tdepBKG(1++)0+_rho_pi_S__b",1.4795);
	Chi2.setParameter("tdepBKG(1++)0+_rho_pi_S__c0",-4.3912);
	Chi2.setParameter("tdepBKG(1++)0+_rho_pi_S__c1",12.770);
	Chi2.setParameter("tdepBKG(1++)0+_rho_pi_S__c2",-23.516);
	Chi2.setParameter("cohBKG(1++)0+_f0(980)_pi_P__b",-5.8146);
	Chi2.setParameter("cohBKG(1++)0+_rho_pi_D__b",2.5961);
	Chi2.setParameter("M_pi(1800)",1.8005);
	Chi2.setParameter("G_pi(1800)",0.21337);
	Chi2.setParameter("cohBKG(0-+)0+_f0(980)_pi_S__b",-3.3995);
	Chi2.setParameter("M_a2(1320)",1.3136);
	Chi2.setParameter("G_a2(1320)",0.11089);
	Chi2.setParameter("M_a2(1700)",1.6704);
	Chi2.setParameter("G_a2(1700)",0.41407);
	Chi2.setParameter("cohBKG(2++)1+_f2_pi_P__b",-0.77708);
	Chi2.setParameter("tdepBKG(2++)1+_rho_pi_D__b",0.66388);
	Chi2.setParameter("tdepBKG(2++)1+_rho_pi_D__c0",-1.6812);
	Chi2.setParameter("tdepBKG(2++)1+_rho_pi_D__c1",0.51046);
	Chi2.setParameter("tdepBKG(2++)1+_rho_pi_D__c2",2.5439);
	Chi2.setParameter("cohBKG(2++)2+_rho_pi_D__b",-0.89233);
	Chi2.setParameter("M_pi2(1670)",1.6521);
	Chi2.setParameter("G_pi2(1670)",0.30451);
	Chi2.setParameter("M_pi2(1880)",1.8335);
	Chi2.setParameter("G_pi2(1880)",0.32157);
	Chi2.setParameter("cohBKG(2-+)0+_f2_pi_D__b",2.4491);
	Chi2.setParameter("tdepBKG(2-+)0+_f2_pi_S__b",2.0204);
	Chi2.setParameter("tdepBKG(2-+)0+_f2_pi_S__c0",-3.7016);
	Chi2.setParameter("tdepBKG(2-+)0+_f2_pi_S__c1",3.1433);
	Chi2.setParameter("tdepBKG(2-+)0+_f2_pi_S__c2",-2.9780);
	Chi2.setParameter("tdepBKG(2-+)0+_rho_pi_F__b",-1.3158);
	Chi2.setParameter("tdepBKG(2-+)0+_rho_pi_F__c0",0.90445);
	Chi2.setParameter("tdepBKG(2-+)0+_rho_pi_F__c1",-0.18399E-01);
	Chi2.setParameter("tdepBKG(2-+)0+_rho_pi_F__c2",2.0252);
	Chi2.setParameter("cohBKG(2-+)1+_f2_pi_S__b",-0.92364);
	Chi2.setParameter("M_a4(2040)",1.9256);
	Chi2.setParameter("G_a4(2040)",0.37209);
	Chi2.setParameter("cohBKG(4++)1+_f2_pi_F__b",-0.18977);
	Chi2.setParameter("cohBKG(4++)1+_rho_pi_G__b",-4.3411);
///////////////////////////////////////////////////////////////////////
	std::cout<< "nPar:"<<nPar<<"; nCpl:"<<nCpl<<"; nBra:"<<nBra<<" => nTot:"<<nTot<<std::endl;

	Chi2.fit();
//	Chi2.conjugate();
//	Chi2.initCouplings();
//	print_vector(Chi2.getParameters());
	double params[] = {1.41204, -1.06975, -0.00748951, 0.0380807, 3.3125, -2.89382, -1.54082, 0.803489, 0.00969979, -0.0418916, -3.49037, 1.92446, -0.365443, -1.64996, 0.0421912, 0.00731507, -0.702389, -3.52486, 1.56651, 0.473571, -0.0379001, 0.0522614, 2.9242, 1.18575, -1.37183, 0.801453, -0.0181871, -0.0704837, -2.56476, 0.971844, 1.52324, -0.234316, -0.0175631, 0.0918839, 2.27319, 0.236614, -1.429, -0.3788, 0.0478653, -0.0855039, -1.44503, -1.0808, 0.280363, 1.33411, -0.0942537, -0.0137978, -0.556254, 1.30749, -0.366226, -1.09952, 0.0941395, -0.01343, 0.7383, -0.913639, -0.0014982, -0.679641, 0.0971191, -0.0280913, 1.31206, 0.661233, -0.265434, 0.216036, -0.0225412, -0.0581884, 1.10534, -0.86583, 1.2794, 0.40958, 2.0108, 0.32212, 1.4795, -4.3912, 12.77, -23.516, 1.4093, 0.13925, -5.8146, 2.5961, 1.8005, 0.21337, -3.3995, 1.3136, 0.11089, 1.6704, 0.41407, -0.77708, 0.66388, -1.6812, 0.51046, 2.5439, -0.89233, 1.6521, 0.30451, 1.8335, 0.32157, 2.4491, 2.0204, -3.7016, 3.1433, -2.978, -1.3158, 0.90445, -0.018399, 2.0252, -0.92364, 1.9256, 0.37209, -0.18977, -4.3411, 1.001, 0.00132763, 0.829259, 0.0181491, -15.3911, 2.69214, 119.216, -83.8281, 1460.37, -348.682, -412.197, 6.58593, 167.65, -7.31843, -7.05273, 0.651271, 1546.29, -2450.71, -2127.02, 3600.22, 87.4355, -98.3871, -43.6377, 7.59307, 310.204, 32.4695, 139.365, -40.1092, 9670.56, -4562.28, 909.894, -355.194};

//	Chi2.open_output("../compare/chi2_f.dat");
//	std::cout<<"chi2: "<<Chi2(params)<<std::endl;
//	Chi2.close_output();
//	for (int i=2*nCpl;i<nTot;i++){
//		std::cout <<"RelPar('"<<Chi2.getParName(i)<<"'):"<<std::endl;
//		Chi2.relPar(i);
//		std::cout<<Chi2.fit()<<std::endl;
//	};
//	print_vector(Chi2.getParameters());
//	Chi2.printStatus();
	return 0;
};

