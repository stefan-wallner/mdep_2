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

#ifdef ADOL_ON
#include "adolc/adolc.h"
#endif//ADOL_ON

#include"matrix_utilities.h"

minimize::minimize(): 
	_method(), 
	_init(false),
	_randRange(100.), 
	_maxFunctionCalls(1000000),
	_maxIterations(100000),
	_tolerance(1.),
	_minStepSize(0.0001){};
#ifdef USE_YAML
//########################################################################################################################################################
///Constructor from YAML file
minimize::minimize(
							std::string 					card):
	_method(card), 
	_init(false),
	_randRange(100.), 
	_maxFunctionCalls(1000000),
	_maxIterations(100000),
	_tolerance(1.),
	_minStepSize(0.0001){

	YAML::Node Ycard   = YAML::LoadFile(card);
	std::cout<<"Load fitter definitions"<<std::endl;
	loadFitterDefinitions(Ycard);
	std::cout<<"Fitter definitions loaded\nFinish setup"<<std::endl;
	finish_setUp();
	std::cout<<"Setup finished"<<std::endl;
};
#endif//USE_YAML
//########################################################################################################################################################
///Actual call for fitter. At the moment the instance is copied and fitted, this might be improved...
double minimize::fit(){ 

	print_vector(_released);
	_f=ROOT::Math::Functor(_method,_method.nTot());
	if(_init){
		_min->Minimize();
		const double *xs = _min->X();
		std::vector<double> best_par(_method.nTot());
		for (int i=0;i<_method.nTot();i++){
			_method.setParameter(i,xs[i]);
		};
		return _method(xs);
	}else{
		std::cerr<<"Error: Fitter not initialized"<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};
//########################################################################################################################################################
///Finds start values for couplings and branchings
void minimize::initCouplings(){ 

	std::cout<<"Initialize couplings"<<std::endl;
	int cpls = _method.nBrCplAnc();
	setRandomCpl(); // Set random couplings
	for(int tbin=0; tbin<_method.Waveset()->nTbin();tbin++){ // Switch off all t' bins
		_method.Waveset()->setEvalTbin(tbin,false);
	};
	_method.setUseBranch(false); // Do not use branchings at first
	for(int tbin=0; tbin<_method.Waveset()->nTbin();tbin++){ // Find cpl for each t' bin
		_method.Waveset()->setEvalTbin(tbin,true);
		std::cout<<"tBin #"<<tbin<<std::endl;
		for (int i =0;i<2*_method.nCpl();i++){
			fixPar(i);
		};
		for (int i=0;i<2*cpls;i++){
			relPar(2*cpls*tbin+i);
		};
		double onetbinchi2 = fit();
		std::cout <<"... Chi2 = "<<onetbinchi2<<std::endl;
		_method.Waveset()->setEvalTbin(tbin,false);
	};
	for(int tbin=0;tbin<_method.Waveset()->nTbin();tbin++){ // Switch on all t' bins
		_method.Waveset()->setEvalTbin(tbin,true);
	};
	std::vector<std::complex<double> > couplings(_method.nCpl());
	std::vector<double> par(_method.nPar());
	for (int i=0;i<_method.nCpl();i++){
		couplings[i] = std::complex<double>(_method.parameters()[2*i],_method.parameters()[2*i+1]);
	};
	for (int i=0;i<_method.nPar();i++){
		par[i] = _method.parameters()[2*_method.nCpl()+i];
	};
	std::vector<double> iso_par(_method.nIso());
	for (int i=0;i<_method.nIso();i++){
		iso_par[i] = _method.parameters()[2*_method.nCpl()+_method.nPar()+2*_method.nBra()+i];
	};
	std::cout << "Total with _method.EvalAutoCpl() (For consistency check): "<< _method.EvalAutoCpl(&couplings[0],&par[0],&iso_par[0])<<std::endl;
	_method.setUseBranch(true);
	if (_method.nBra()>0){
		std::vector<std::complex<double> > bra = _method.get_branchings(couplings,par,iso_par);
		// branchCouplingsToOne(); // Set all coupled couplings to one, since all should be in the branchings right now // Somehow Chi2 is better, when this is not done
		for (unsigned int i=0;i<bra.size();i++){ // Set found branchings
			setParameter(2*_method.nCpl()+_method.nPar()+2*i ,bra[i].real());
			setParameter(2*_method.nCpl()+_method.nPar()+2*i+1,bra[i].imag());
		};
		std::cout << "With the found branchings, Chi2(...)="<< _method.EvalAutoCplBranch(&bra[0],&couplings[0],&par[0],&iso_par[0])<<" ('_method.EvalAutoCplBranch(...)')"<<std::endl; //[0]//
		for (int i =0;i<2*_method.nCpl();i++){ // Fix couplings
			fixPar(i);
		};
		for (int i=0;i<2*_method.nBra();i++){ // Rel Branchings
			relPar(2*_method.nCpl()+_method.nPar()+i);
		};
		fit();
		std::cout<<"Couplings and branchings"<<std::endl;
		for (int i =0;i<2*_method.nCpl();i++){ // Rel couplings
			relPar(i);
		};
		fit();
	}else{
		for(int i=0;i<_method.nCpl();i++){
			relPar(2*i);
			relPar(2*i+1);
		};
	};
	std::cout<<"Total: "<<_method(_min->X())<<std::endl;
	std::cout<<"Couplings and branchings found"<<std::endl;
	std::cout<<"Setting automatic limits for couplings and branchings"<<std::endl;
	for (int i=0;i<_method.nCpl();i++){
		double val = std::max(_method.parameters()[2*i]*_method.parameters()[2*i],_method.parameters()[2*i+1]*_method.parameters()[2*i+1]);
		val = pow(val,.5);
		_method.setParLimits(2*i  ,3*val,-3*val);
		_method.setParLimits(2*i+1,3*val,-3*val);
	};
	int par_bef = 2*_method.nCpl() +_method.nPar();
	for (int i=0;i<_method.nBra();i++){
		double val = std::max(_method.parameters()[par_bef+2*i]*_method.parameters()[par_bef+2*i],_method.parameters()[par_bef+2*i+1]*_method.parameters()[par_bef+2*i+1]);
		val = pow(val,.5);
		_method.setParLimits(par_bef+2*i  ,3*val,-3*val);
		_method.setParLimits(par_bef+2*i+1,3*val,-3*val);
	};
};
//########################################################################################################################################################
///Call lower 'setParamter' and sets internal definitions
void minimize::setParameter(
							int 						i, 	// # of parameter
							double 						par){

	_method.setParameter(i,par);
	reload_par_definitions(i);
};
//########################################################################################################################################################
///Sets a parameter by name
void minimize::setParameter(
							std::string 					name, 
							double 						par){

	int number = _method.getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_method.nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		setParameter(number,par);
	};
};
//########################################################################################################################################################
///Sets all parameters at once
void minimize::setParameters(
							std::vector<double> 				pars){

	_method.setParameters(pars);
	reload_par_definitions();
};
//########################################################################################################################################################
///Sets parameter step size
void minimize::setStepSize(
							int 						i, 	// # of parameter
							double 						step){

	if (_step_sizes.size() < _method.nTot()){
		_step_sizes = std::vector<double>(_method.nTot(),_minStepSize);
		reload_par_definitions();
	};
	_step_sizes[i] = step;
	reload_par_definitions(i);
};

//########################################################################################################################################################
///Sets the step size by name
void minimize::setStepSize(
							std::string 					name, 
							double 						par){

	int number = _method.getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_method.nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		setStepSize(number,par);
	};
};
//########################################################################################################################################################
///Sets all step sizes
void minimize::setStepSizes(
							std::vector<double> 				steps){

	if (steps.size() >= _method.nTot()){
		_step_sizes=steps;
		reload_par_definitions();
	}else{
		std::cerr<<"Error: Input steps.size() too small"<<std::endl;
	};
};
//########################################################################################################################################################
///Sets the range for generating random couplings and branchings
void minimize::setRandRange(
							double 						range){

	_randRange=range;
};
//########################################################################################################################################################
///Release parameter by number
void minimize::relPar(
							int 						i){	// # of parameter

	if (_released.size() == _method.nTot()){
		_released[i] = true;
		reload_par_definitions(i);
	}else{
		std::cout<<_method.nTot()<<" "<<_method.parameters().size()<<" "<<_released.size()<<std::endl;
		std::cerr<<"Error: _released is not initialized."<<std::endl;
	};
};
//########################################################################################################################################################
///Fix parameter by number
void minimize::fixPar(
							int 						i){	// # of parameter

	if (_released.size() == _method.nTot()){
		_released[i] = false;
		reload_par_definitions(i);
	}else{
		std::cout<<_method.nTot()<<" "<<_method.parameters().size()<<" "<<_released.size()<<std::endl;
		std::cerr<<"Error: _released is not initialized."<<std::endl;
	};
};
//########################################################################################################################################################
///Release parameter by name
void minimize::relPar(
							std::string 					name){

	int number = _method.getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_method.nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		relPar(number);
	};
};
//########################################################################################################################################################
///Fix parameter by name
void minimize::fixPar(
							std::string 					name){

	int number = _method.getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if (_method.nTot() <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
	}else{
		fixPar(number);
	};
};
//########################################################################################################################################################
///Prints the internal status
void minimize::printStatus(){ 

	_method.printStatus();
	std::cout<<std::endl<<"_best_par"<<std::endl;
	print_vector(_best_par);
	std::cout<<std::endl<<"_randRange: "<<_randRange<<std::endl;
	std::cout<<std::endl<<"_minStepSize: "<<_minStepSize<<std::endl;
	std::cout<<std::endl<<"_init: "<<_init<<std::endl;
	std::cout<<std::endl<<"_maxFunctionCalls: "<<_maxFunctionCalls<<std::endl;
	std::cout<<std::endl<<"_maxIterations: "<<_maxIterations<<std::endl;
	std::cout<<std::endl<<"_tolerance: "<<_tolerance<<std::endl;
	std::cout<<std::endl<<"_step_sizes"<<std::endl;
	print_vector(_step_sizes);
	std::cout<<std::endl<<"_released:"<<std::endl;
	print_vector(_released);
};
//########################################################################################################################################################
///Updates internal definitions
void minimize::update_definitions(){ 
	_method.update_definitions();
	std::vector<bool> rels;
	for (int i=0;i<_method.nCpl();i++){
		rels.push_back(true);
		rels.push_back(true);
	};
	for (int i=0;i<_method.nPar();i++){
		rels.push_back(false);
	};
	for (int i=0;i<_method.nBra();i++){
		rels.push_back(false);
		rels.push_back(false);
	};
	for (int i=0;i<_method.nIso();i++){
		rels.push_back(false);
	};
	_released = rels;
	if(_init){
		_min->SetTolerance(_tolerance);
		_min->SetMaxFunctionCalls(_maxFunctionCalls);
		_min->SetMaxIterations(_maxIterations);
	};
};
//########################################################################################################################################################
///Update the parameter definitions for parameter number mara_peter, if mara_peter ==-1, then all are updated
void minimize::reload_par_definitions(
							int 						mara_peter){ 

	int uLim = 0;
	int oLim = _method.nTot();
	if (mara_peter > -1){
		uLim = mara_peter;
		oLim = mara_peter+1;
	};
	if((*_method.lower_parameter_limits()).size() != _method.parameters().size()){
		_method.init_lower_limits(_method.parameters().size());
	};
	if((*_method.upper_parameter_limits()).size() != _method.parameters().size()){
		_method.init_upper_limits(_method.parameters().size());
	};
	if(_init){
		for(int i=uLim;i<oLim;i++){
			if(_released[i]){
				if((*_method.lower_parameter_limits())[i] < (*_method.upper_parameter_limits())[i]){
					_min->SetLimitedVariable(i,(*_method.parNames())[i],_method.parameters()[i],_step_sizes[i],(*_method.lower_parameter_limits())[i],(*_method.upper_parameter_limits())[i]);
				}else{
					_min->SetVariable(i,(*_method.parNames())[i],_method.parameters()[i],_step_sizes[i]);
				};
			}else{
				_min->SetFixedVariable(i,(*_method.parNames())[i],_method.parameters()[i]);
			};
		};
	};
};
//########################################################################################################################################################
///Initializes the fitter
bool minimize::initialize(std::string s1, std::string s2){ 

	if (_method.parameters().size()<_method.nTot()){
		std::cout<<"_method.parameters().size() < _method.nTot(). Abort initialize(initialization."<<std::endl;
		return false;
	};
	if ((*_method.parNames()).size()<_method.nTot()){
		std::cout<<"(*_method.parNames()).size() < _method.nTot(). Abort initialization."<<std::endl;
		return false;
	};
	if (_step_sizes.size()<_method.nTot()){
		std::cout<<"_step_sizes.size() < _method.nTot(). Abort initialization."<<std::endl;
		return false;
	};
	if (_released.size()<_method.nTot()){
		std::cout<<"_released.size() < _method.nTot(). Abort initialization."<<std::endl;
		return false;
	};
	_min = ROOT::Math::Factory::CreateMinimizer(s1,s2);
	_min->SetMaxFunctionCalls(_maxFunctionCalls);
	_min->SetMaxIterations(_maxIterations);
	_min->SetTolerance(_tolerance);
	_f=ROOT::Math::Functor(_method,_method.nTot());
	_min->SetFunction(_f);
	_init = true;
	update_definitions();
	reload_par_definitions();
	return true;
};
//########################################################################################################################################################
///Randomize couplings within _randRange
void minimize::setRandomCpl(){

	for(int i=0; i<2*_method.nCpl();i++){
		if (_released[i]){
			_method.setParameter(i,2*((double)rand()/RAND_MAX-0.5)*_randRange);
		};
	};
	reload_par_definitions();
};
//########################################################################################################################################################
///Randomize branchings within _randRange
void minimize::setRandomBra(){ 

	for (int i=0;i<2*_method.nBra();i++){
		if (_released[2*_method.nCpl()+_method.nPar()+i]){
			_method.setParameter(2*_method.nCpl()+_method.nPar()+i, 2*((double)rand()/RAND_MAX-0.5)*_randRange);
		};
	};
	reload_par_definitions();
};
//########################################################################################################################################################
///Sets some internal vaiables accordingly
void minimize::finish_setUp(){

	std::vector<int> first_branch = _method.Waveset()->getFirstBranch();
	if (_released.size() != _method.parameters().size()){
	//		std::cout<<"Warning: No paramter status set, releasing couplings, fixing all others."<<std::endl;
		std::vector<bool> std_rel(_method.parameters().size(),false);
		for (int i=0;i<2*_method.nCpl();i++){
			std_rel[i]=true;
		};
		_released = std_rel;
	};
	for (int i=0; i<first_branch.size();i++){
		setParameter(2*_method.nCpl()+_method.nPar()+2*first_branch[i]  ,1.); // Re(Br)
		setParameter(2*_method.nCpl()+_method.nPar()+2*first_branch[i]+1,0.); // Im(Br)
		fixPar(2*_method.nCpl()+_method.nPar()+2*first_branch[i]  );
		fixPar(2*_method.nCpl()+_method.nPar()+2*first_branch[i]+1);
	};
	setRandomCpl();
	setRandomBra();
	initialize();
};
#ifdef USE_YAML
//########################################################################################################################################################
///Loads some fitter definitions from a YAML file
void minimize::loadFitterDefinitions(
							YAML::Node 					&waveset){

	if (waveset["min_step_size"]){
		_minStepSize = waveset["min_step_size"].as<double>();
	};
	int nPar = _method.Waveset()->getNpar();
	int nCpl = _method.getNanc();
	for (int par = 0;par<nPar;par++){
		setStepSize(2*nCpl+par,std::max(_minStepSize, 0.0001*fabs(_method.parameters()[2*nCpl+par])));
	};
	if (waveset["real_anchor_cpl"]){ /// Also include this here, so no extra method is necessary
		if(waveset["real_anchor_cpl"].as<bool>()){
			setParameter(1,0.);
			fixPar(1);
		};
	};
};
#endif//USE_YAML
//########################################################################################################################################################
///Cube required by the MultiNest package
void minimize::cube(					double						*in)						const{
	
	for (int i=0;i<_method.nTot();i++){
		in[i] = (1-in[i])*(*_method.lower_parameter_limits())[i]+ in[i]*(*_method.upper_parameter_limits())[i];
	};
};
//########################################################################################################################################################
