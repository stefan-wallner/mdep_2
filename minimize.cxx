#include "minimize.h"

#include<vector>
#include<complex>
#include<string>
#include<iostream>
#include<fstream>
#include<cstdlib>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "adolc/adolc.h"

#include"invert33.h"

minimize::minimize(): anchor_t(), _init(false),_randRange(100.), _maxFunctionCalls(1000000),_maxIterations(100000),_tolerance(1.),_minStepSize(0.0001){};


#ifdef USE_YAML
///Constructror from YAML files
minimize::minimize(std::string card, std::string waves, std::string parametrizations):
	anchor_t(card,waves,parametrizations), _init(false),_randRange(100.), _maxFunctionCalls(1000000),_maxIterations(100000),_tolerance(1.),_minStepSize(0.0001){
	YAML::Node Ycard   = YAML::LoadFile(card);
	loadFitterDefinitions(Ycard);
	finish_setUp();
};
#endif//USE_YAML

std::vector<bool> minimize::getReleased(){// Self explanatory
	return _released;
};


void minimize::finish_setUp(){
	std::vector<int> first_branch = getFirstBranch();
	if (_released.size() != _parameters.size()){
		std::cout<<"Warning: No paramter status set, releasing couplings, fixing all others."<<std::endl;
		std::vector<bool> std_rel(_parameters.size(),false);
		for (int i=0;i<2*_nCpl;i++){
			std_rel[i]=true;
		};
		_released = std_rel;
	};
	for (int i=0; i<first_branch.size();i++){
		setParameter(2*_nCpl+_nPar+2*first_branch[i]  ,1.); // Re(Br)
		setParameter(2*_nCpl+_nPar+2*first_branch[i]+1,0.); // Im(Br)
		fixPar(2*_nCpl+_nPar+2*first_branch[i]  );
		fixPar(2*_nCpl+_nPar+2*first_branch[i]+1);
	};
	setRandomCpl();
	setRandomBra();
	initialize();
};


void minimize::setStepSize(std::string name, double par){// Self explanatory
	int number = getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_nTot <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		setStepSize(number,par);
	};
};

void minimize::fixPar(std::string name){// Self explanatory
	int number = getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_nTot <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		fixPar(number);
	};
};

void minimize::relPar(std::string name){// Self explanatory
	int number = getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_nTot <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		relPar(number);
	};
};






void minimize::setStepSize(int i, double step){
	if (_step_sizes.size() < _nTot){
		_step_sizes = std::vector<double>(_nTot,_minStepSize);
		reload_par_definitions();
	};
	_step_sizes[i] = step;
	reload_par_definitions(i);
};

void minimize::setStepSizes(std::vector<double> steps){ // Self explanatory
	if (steps.size() >= _nTot){
		_step_sizes=steps;
		reload_par_definitions();
	}else{
		std::cerr<<"Error: Input steps.size() too small"<<std::endl;
	};
};

void minimize::fixPar(int i){ // fix parameter by number
	if (_released.size() == _nTot){
		_released[i] = false;
		reload_par_definitions(i);
	}else{
		std::cout<<_nTot<<" "<<_parameters.size()<<" "<<_released.size()<<std::endl;
		std::cerr<<"Error: _released is not initialized."<<std::endl;
	};
};

void minimize::relPar(int i){ // Release parameter by number
	if (_released.size() == _nTot){
		_released[i] = true;
		reload_par_definitions(i);
	}else{
		std::cout<<_nTot<<" "<<_parameters.size()<<" "<<_released.size()<<std::endl;
		std::cerr<<"Error: _released is not initialized."<<std::endl;
	};
};

void minimize::reload_par_definitions(int mara_peter){ // Update the parameter definitions for parameter number mara_peter, if mara_peter ==-1, then all are updated
	int uLim = 0;
	int oLim = _nTot;
	if (mara_peter > -1){
		uLim = mara_peter;
		oLim = mara_peter+1;
	};
	if(_lower_parameter_limits.size() != _parameters.size()){
		_lower_parameter_limits=std::vector<double>(_parameters.size(),1.);
	};
	if(_upper_parameter_limits.size() != _parameters.size()){
		_upper_parameter_limits=std::vector<double>(_parameters.size(),0.);
	};
	if(_init){
		for(int i=uLim;i<oLim;i++){
			if(_released[i]){
				if(_lower_parameter_limits[i] < _upper_parameter_limits[i]){
					_min->SetLimitedVariable(i,_parNames[i],_parameters[i],_step_sizes[i],_lower_parameter_limits[i],_upper_parameter_limits[i]);
				}else{
					_min->SetVariable(i,_parNames[i],_parameters[i],_step_sizes[i]);
				};
			}else{
				_min->SetFixedVariable(i,_parNames[i],_parameters[i]);
			};
		};
	};
};


std::string minimize::className(){
	return "minimize";
};


bool minimize::initialize(std::string s1, std::string s2){ // Initializes the fitter (Has to be done at least once)
	if (_parameters.size()<_nTot){
		std::cout<<"_parameters.size() < _nTot. Abort initialize(initialization."<<std::endl;
		return false;
	};
	if (_parNames.size()<_nTot){
		std::cout<<"_parNames.size() < _nTot. Abort initialization."<<std::endl;
		return false;
	};
	if (_step_sizes.size()<_nTot){
		std::cout<<"_step_sizes.size() < _nTot. Abort initialization."<<std::endl;
		return false;
	};
	if (_released.size()<_nTot){
		std::cout<<"_released.size() < _nTot. Abort initialization."<<std::endl;
		return false;
	};
	_min = ROOT::Math::Factory::CreateMinimizer(s1,s2);
	_min->SetMaxFunctionCalls(_maxFunctionCalls);
	_min->SetMaxIterations(_maxIterations);
	_min->SetTolerance(_tolerance);
	// _f=ROOT::Math::Functor(fcn_pointer,_nTot);
	_f=ROOT::Math::Functor((*this),_nTot);
	_min->SetFunction(_f);
	_init = true;
	update_definitions();
	reload_par_definitions();
	return true;
};

void minimize::update_definitions(){ // Updates definitions
	anchor_t::update_definitions();
	std::vector<bool> rels;
	for (int i=0;i<_nCpl;i++){
		rels.push_back(true);
		rels.push_back(true);
	};
	for (int i=0;i<_nPar;i++){
		rels.push_back(false);
	};
	for (int i=0;i<_nBra;i++){
		rels.push_back(false);
		rels.push_back(false);
	};
	_released = rels;
	if(_init){
		_min->SetTolerance(_tolerance);
		_min->SetMaxFunctionCalls(_maxFunctionCalls);
		_min->SetMaxIterations(_maxIterations);
	};
};

double minimize::fit(){ // Actual call for fitter. At the moment the instance is copied and fitted, this might be improved...
	print_vector(_released);
	_f=ROOT::Math::Functor((*this),_nTot);
	if(_init){
		_min->Minimize();
		const double *xs = _min->X();
		std::vector<double> best_par(_nTot);
		for (int i=0;i<_nTot;i++){
			_parameters[i] = xs[i];
		};
		return (*this)(xs);
	}else{
		std::cerr<<"Error: Fitter not initialized"<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};

void minimize::printStatus(){ // Prints the internal status
	anchor_t::printStatus();
	std::cout<<std::endl<<"_init: "<<_init<<std::endl;
	std::cout<<std::endl<<"_maxFunctionCalls: "<<_maxFunctionCalls<<std::endl;
	std::cout<<std::endl<<"_maxIterations: "<<_maxIterations<<std::endl;
	std::cout<<std::endl<<"_tolerance: "<<_tolerance<<std::endl;
	std::cout<<std::endl<<"_randRange: "<<_randRange<<std::endl;
	std::cout<<std::endl<<"_best_par"<<std::endl;
	print_vector(_best_par);
	std::cout<<std::endl<<"_minStepSize: "<<_minStepSize<<std::endl;
	std::cout<<std::endl<<"_step_sizes"<<std::endl;
	print_vector(_step_sizes);
	std::cout<<std::endl<<"_released:"<<std::endl;
	print_vector(_released);
};

void minimize::setRandRange(double range){
	_randRange=range;
};

void minimize::setRandomCpl(){// Randomize couplings within _randRange
	for(int i=0; i<2*_nCpl;i++){
		if (_released[i]){
			_parameters[i] = 2*((double)rand()/RAND_MAX-0.5)*_randRange;
		};
	};
	reload_par_definitions();
};

void minimize::setRandomBra(){ // Randomize branchings within _randRange
	for (int i=0;i<2*_nBra;i++){
		if (_released[2*_nCpl+_nPar+i]){
			_parameters[2*_nCpl+_nPar+i] = 2*((double)rand()/RAND_MAX-0.5)*_randRange;
		};
	};
	reload_par_definitions();
};

void minimize::initCouplings(){ // Finds 'good' values for couplings and branchings
	std::cout<<"Initialize couplings"<<std::endl;
	int cpls = _nBrCplAnc;
	setRandomCpl(); // Set random couplings
	for(int tbin=0; tbin<_nTbin;tbin++){ // Switch off all t' bins
		setEvalTbin(tbin,false);
	};
	_useBranch = false; // Do not use branchings at first
	for(int tbin=0; tbin<_nTbin;tbin++){ // Find cpl for each t' bin
		setEvalTbin(tbin,true);
		std::cout<<"tBin #"<<tbin<<std::endl;
		for (int i =0;i<2*_nCpl;i++){
			fixPar(i);
		};
		for (int i=0;i<2*cpls;i++){
			relPar(2*cpls*tbin+i);
		};
		double onetbinchi2 = fit();
		std::cout <<"... Chi2 = "<<onetbinchi2<<std::endl;
		setEvalTbin(tbin,false);
	};
	for(int tbin=0;tbin<_nTbin;tbin++){ // Switch on all t' bins
		setEvalTbin(tbin,true);
	};
	std::vector<std::complex<double> > couplings(_nCpl);
	std::vector<double> par(_nPar);
	for (int i=0;i<_nCpl;i++){
		couplings[i] = std::complex<double>(_parameters[2*i],_parameters[2*i+1]);
	};
	for (int i=0;i<_nPar;i++){
		par[i] = _parameters[2*_nCpl+i];
	};
	std::vector<double> iso_par;//herre
	std::cout << "Total with EvalAutoCpl() (For consistency check): "<< EvalAutoCpl(couplings,par,iso_par)<<std::endl;
	_useBranch=true;
	if (_nBra>0){
		std::vector<double> iso_par;//herre
		std::vector<std::complex<double> > bra = get_branchings(couplings,par,iso_par);
		// branchCouplingsToOne(); // Set all coubled couplings to one, since all should be in the branchings right now // Somehow Chi2 is better, when this is not done
		for (unsigned int i=0;i<bra.size();i++){ // Set found branchings
			setParameter(2*_nCpl+_nPar+2*i ,bra[i].real());
			setParameter(2*_nCpl+_nPar+2*i+1,bra[i].imag());
		};
		std::cout << "With the found branchings, Chi2(...)="<< EvalAutoCplBranch(bra,couplings,par,iso_par)<<" ('EvalAutoCplBranch(...)')"<<std::endl;
		for (int i =0;i<2*_nCpl;i++){ // Fix couplings
			fixPar(i);
		};
		for (int i=0;i<2*_nBra;i++){ // Rel Branchings
			relPar(2*_nCpl+_nPar+i);
		};
		fit();
		std::cout<<"Couplings and branchings"<<std::endl;
		for (int i =0;i<2*_nCpl;i++){ // Rel couplings
			relPar(i);
		};
		fit();
	}else{
		for(int i=0;i<_nCpl;i++){
			relPar(2*i);
			relPar(2*i+1);
		};
	};
	std::cout<<"Total: "<<(*this)(_min->X())<<std::endl;
	std::cout<<"Couplings and branchings found"<<std::endl;
};






void minimize::setParameter(int i, double par){
	anchor_t::setParameter(i,par);
	reload_par_definitions(i);
};
void minimize::setParameter(std::string name, double par){
	int number = getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_nTot <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		setParameter(number,par);
	};
};

void minimize::setParameters(std::vector<double> pars){
	anchor_t::setParameters(pars);
	reload_par_definitions();
};
#ifdef USE_YAML
void minimize::loadFitterDefinitions(YAML::Node &waveset){
	if (waveset["min_step_size"]){
		_minStepSize = waveset["min_step_size"].as<double>();
	};
	int nPar = getNpar();
	int nCpl = getNanc();
	for (int par = 0;par<nPar;par++){
		setStepSize(2*nCpl+par,max(_minStepSize, 0.0001*fabs(_parameters[2*nCpl+par])));
	};
	if (waveset["real_anchor_cpl"]){ /// Also include this here, so no extra method is necessary
		if(waveset["real_anchor_cpl"].as<bool>()){
			setParameter(1,0.);
			fixPar(1);
		};
	};
};

#endif//USE_YAML

