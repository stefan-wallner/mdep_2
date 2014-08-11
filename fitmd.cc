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

	std::cout<< "nPar:"<<nPar<<"; nCpl:"<<nCpl<<"; nBra:"<<nBra<<" => nTot"<<nTot<<std::endl;
	Chi2.initCouplings();
	for (int i=2*nCpl;i<nTot;i++){
		std::cout <<"RelPar('"<<Chi2.getParName(i)<<"'):"<<std::endl;
		Chi2.relPar(i);
		std::cout<<Chi2.fit()<<std::endl;
	};
	print_vector(Chi2.getParameters());
	return 0;
};

