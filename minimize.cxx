#include "minimize.h"
#include<vector>
#include<complex>
#include<string>
#include<iostream>
#include<fstream>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include"invert33.h"

minimize::minimize(): anchor_t(), _init(false), _nOut(1000), _count(0), _maxFunctionCalls(1000000),_maxIterations(100000),_tolerance(1.),_minStepSize(0.0001),_randRange(100.),_useBranch(true){};

std::vector<std::string> minimize::getParNames(){// Self explanatory
	return _parNames;
};

std::vector<bool> minimize::getReleased(){// Self explanatory
	return _released;
};

std::vector<double> minimize::getParameters(){// Self explanatory
	return _parameters;
};

void minimize::setParameter(std::string name, double par){// Self explanatory
	int number = getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_nTot <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		setParameter(number,par);
	};
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

int minimize::getParNumber(std::string name){ // Gets the number for a given name, if the name does not exist, return -1
	for (int i=0;i<_parNames.size();i++){
		if (_parNames[i] == name){
			return i;
		};
	};
	return -1;
};

void minimize::setParameter(int i, double par){
	if(_parameters.size() < _nTot){
		_parameters = std::vector<double>(_nTot,0.);
		reload_par_definitions();
	};
	_parameters[i] = par;
	reload_par_definitions(i);
};

void minimize::setParameters(std::vector<double> pars){
	if(pars.size()>=_nTot){
		_parameters = pars;
		reload_par_definitions();
	}else{
		std::cerr<<"Error: Input pars.size() too small"<<std::endl;
	};
};

void minimize::setParLimits(int par, double upper, double lower){
	if(_lower_parameter_limits.size() != _parameters.size()){
		_lower_parameter_limits=std::vector<double>(_parameters.size(),std::numeric_limits<double>::quiet_NaN());
	};
	if(_upper_parameter_limits.size() != _parameters.size()){
		_upper_parameter_limits=std::vector<double>(_parameters.size(),std::numeric_limits<double>::quiet_NaN());
	};
	_upper_parameter_limits[par]=upper;
	_lower_parameter_limits[par]=lower;
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
		std::cerr<<"Error: _released is not initialized."<<std::endl;
	};
};

void minimize::relPar(int i){ // Release parameter by number
	if (_released.size() == _nTot){
		_released[i] = true;
		reload_par_definitions(i);
	}else{
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
		_lower_parameter_limits=std::vector<double>(_parameters.size(),std::numeric_limits<double>::quiet_NaN());
	};
	if(_upper_parameter_limits.size() != _parameters.size()){
		_upper_parameter_limits=std::vector<double>(_parameters.size(),std::numeric_limits<double>::quiet_NaN());
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

double minimize::operator()(){ // Evaluate Chi2 with the _parameters (Call the oter operator)
	double pars[2*_nCpl+_nPar+2*_nBra];
	for (int i=0;i<2*_nCpl+_nPar+2*_nBra;i++){
		pars[i]=_parameters[i];
	};
	return (*this)(pars);
};

double minimize::operator()(const double* xx){ // Call of the operator
	std::vector<std::complex<double> > cpl(_nCpl);
	std::vector<double> par(_nPar);
	std::vector<std::complex<double> > bra(_nBra);
	int count_xx =0;
	for (int i=0;i<_nCpl;i++){ // Build cpl, par, bra from xx
		cpl[i] = std::complex<double>(xx[count_xx],xx[count_xx+1]);
		count_xx+=2;
	};
	for (int i=0;i<_nPar;i++){
		par[i] = xx[count_xx];
		count_xx++;
	};
	for (int i=0;i<_nBra;i++){
		bra[i] = std::complex<double>(xx[count_xx],xx[count_xx+1]);
		count_xx+=2;
	};
	// std::cout<<par[0]<<std::endl;
	double chi2;
	if (_useBranch){ // Evaluate
		chi2 = EvalAutoCplBranch(bra,cpl,par);
	}else{
		chi2 = EvalAutoCpl(cpl,par); // This only works, because of the Automatic coupling finding algorithm, switching off the branchings does not give adiitional couplings
	};
	_count++;
	/* if (chi2 != chi2){ // Check for NaN
	std::cout<<"NaN paramters are:"<<std::endl;
	std::cout<<"bra:"<<std::endl;
	print_vector(bra);
	std::cout<<"cpl:"<<std::endl;
	print_vector(cpl);
	std::cout<<"par:"<<std::endl;
	print_vector(par);
	std::cout << xx << std::endl;
	throw;
	};*/
	if (0==_count%_nOut){ // Write every _nOut evaluation
		std::cout<<"#"<<_count<<": "<<chi2<<std::endl;
	};
	return chi2;
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
	_nTot = getNtotAnc();
	_nPar = getNpar();
	_nCpl = getNanc();
	_nBra = getNbra();
	std::vector<std::string> names;
	std::vector<bool> rels;
	std::stringstream str_count;	
	for (int i=0;i<_nCpl;i++){
		str_count<<i;
		names.push_back("reC"+str_count.str());
		names.push_back("imC"+str_count.str());
		str_count.str("");
		rels.push_back(true);
		rels.push_back(true);
	};	
	for (int i=0;i<_nPar;i++){
		names.push_back(getParameterName(i));
		rels.push_back(false);
	};
	for (int i=0;i<_nBra;i++){
		str_count<<i;
		names.push_back("reB"+str_count.str());
		names.push_back("imB"+str_count.str());
		str_count.str("");
		rels.push_back(false);
		rels.push_back(false);
	};
	_parNames = names;
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
	std::cout<<std::endl<<"_parNames:"<<std::endl;
	print_vector(_parNames);
	std::cout<<std::endl<<"_parameters:"<<std::endl;
	print_vector(_parameters);
	std::cout<<std::endl<<"_released:"<<std::endl;
	print_vector(_released);
	std::cout<<std::endl<<"_upper_parameter_limits:"<<std::endl;
	print_vector(_upper_parameter_limits);
	std::cout<<std::endl<<"_lower_parameter_limits:"<<std::endl;
	print_vector(_lower_parameter_limits);	
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
	std::cout << "Total with EvalAutoCpl() (For consistency check): "<< EvalAutoCpl(couplings,par)<<std::endl;
	_useBranch=true;
	if (_nBra>0){
		std::vector<std::complex<double> > bra = get_branchings(couplings,par);
		// branchCouplingsToOne(); // Set all coubled couplings to one, since all should be in the branchings right now // Somehow Chi2 is better, when this is not done
		for (unsigned int i=0;i<bra.size();i++){ // Set found branchings
			setParameter(2*_nCpl+_nPar+2*i ,bra[i].real());
			setParameter(2*_nCpl+_nPar+2*i+1,bra[i].imag());
		};
		std::cout << "With the found branchings, Chi2(...)="<< EvalAutoCplBranch(bra,couplings,par)<<" ('EvalAutoCplBranch(...)')"<<std::endl;
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

void minimize::writePlots(int tbin,std::string filename){ // Does not work with branchings at the moment
	std::vector<std::complex<double> > cpl(_nCpl);
	std::vector<double> par(_nPar);
	int count_xx =0;
	for (int i=0;i<_nCpl;i++){
		cpl[i] = std::complex<double>(_parameters[count_xx],_parameters[count_xx+1]);
		count_xx+=2;
	};
	for (int i=0;i<_nPar;i++){
		par[i] = _parameters[count_xx];
		count_xx++;
	};
	std::vector<std::vector<double> >plots = getPlots(tbin,cpl,par);
	std::ofstream outfile;
	outfile.open(filename.c_str());
	for (int i=0;i<plots[0].size()/3;i++){
		outfile<<(_binning[i]+_binning[i+1])/2<<" ";
		for (int j=0;j<plots.size();j++){
			outfile<<plots[j][3*i]<<" "<<plots[j][3*i+1]<<" "<<pow(plots[j][3*i+2],-.5)<<" ";
		};
		outfile<<"\n";
	};
	outfile.close();
};

std::string minimize::getParName(int i){ // Self explanatory
	return _parNames[i];
};

void minimize::branchCouplingsToOne(){ // Set the anchor coupligs to one, that are coupled to another function
	for (unsigned int i=0;i<_coupled.size();i++){
		if(_coupled[i] >=0){
			int nCpl = _n_cpls[i];
			if (nCpl <= _nBrCplAnc){
				for (int tbin =0;tbin<_nTbin;tbin++){
					_parameters[2*_nBrCplAnc*tbin+2*nCpl ]=1.;
					_parameters[2*_nBrCplAnc*tbin+2*nCpl+1]=0.;
				};
			};
		};
	};
};
