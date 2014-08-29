#include<iostream>
#include"chi2.h"
#include"../chi_squared/breitWigners.h"
#include"../chi_squared/phaseSpace.h"

#include"invert33.h"

double UPPER_MASS_LIMIT = 2.5;
double LOWER_MASS_LIMIT = 0.5;
int DEFAULT_L =0;

std::vector<std::complex<double> > chi2::amps_class(double m, std::vector<double> &par){ // Give the amplitudes with the 'old' parameter-vector format (couplings anf shape parameters mixed)
	std::vector<double> params = std::vector<double>(_nPar);
	std::vector<std::complex<double> > cpl = std::vector<std::complex<double> >(_nFtw);
	int i_max = _interface.size();
	int npar = 0;
	int ncpl = 0;
	for (unsigned int i=0; i<i_max;i++){
		int mode  = _interface[i];
		if(mode==3){
			params[npar]=par[i];
			npar+=1;
		}else if (mode == 2 or mode == 4){ // Imag parts handeled by 'mode==1', consts handeled by internal constants.
			continue;
		}else if (mode == 1){
			cpl[ncpl]=std::complex<double>(par[i],par[i+1]);
			ncpl+=1;
		};
	};
	return amps(m,cpl,params);
};



chi2::chi2():_nWaves(0),_nFuncs(0),_globalPs(0),_maxNpars(0),_nPar(0),_nFtw(0),_write_out(false){};

template<typename xdouble>
std::vector<std::complex<xdouble> > chi2::amps(double m,std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par){ // Builds the amplitudes for each wave from functions with shape parameters par and cpl at m3pi = m 
	std::vector<std::complex<xdouble> > funcEval = funcs(m,par); // Evalulated BW-functions
	std::vector<std::complex<xdouble> > ampl = std::vector<std::complex<xdouble> >(_nWaves); // Vector of final amplitudes
	std::vector<double> ps = phase_space(m);
	std::complex<xdouble> amp; // Actual wave amplitude
	int upBor=0; // Upper limit for function number
	int loBor=0; // Lower limit for function number
	for (int wave =0;wave<_nWaves;wave++){
		amp = std::complex<xdouble>(0.,0.);
		loBor = upBor;
		upBor = _borders_waves[wave];	
		if (m >= _lowerLims[wave] and m<_upperLims[wave]){ 	
			for (int nFunc = loBor; nFunc<upBor; nFunc++){
				int func = _funcs_to_waves[nFunc]; // Number of function contributing to the actual wave
				amp+=std::complex<xdouble>(ps[wave],0.)*cpl[nFunc]*funcEval[func];
			};
		};
		ampl[wave]=amp;
	};
	return ampl;
};
template std::vector<std::complex<double> > chi2::amps(double m, std::vector<std::complex<double> > &cpl,std::vector<double> &par);

template<typename xdouble>
std::vector<std::complex<xdouble> > chi2::funcs(double m,std::vector<xdouble> &par){ // returns the function values at m3pi = m and shape parameters par
	std::vector<std::complex<xdouble> > f = std::vector<std::complex<xdouble> >(_nFuncs);
	int upPar=0;
	int loPar=0;
	int upConst=0;
	int loConst=0;
	std::vector<xdouble> act_par = std::vector<xdouble>(_maxNpars);
	for (int func=0; func<_nFuncs; func++){
		loPar = upPar;
		upPar = _borders_par[func];
		loConst = upConst;
		upConst = _borders_const[func];
		int pos=0;
		if (m >= _funcLowerLims[func] and m < _funcUpperLims[func]){ // Only calculate needed functions
			int nFunc = _funcs[func];
			for (int i = loPar;i<upPar;i++){
				act_par[pos]=par[i];
				pos+=1;
			};
			for (int i = loConst; i<upConst;i++){
				act_par[pos]=_const[i];
				pos+=1;
			};
			f[func]=bw(m,act_par,nFunc,_L_func[func]);
		}else{
			f[func]=std::complex<xdouble>(0.,0.);
		};
			
	};
	return f;
};
template std::vector<std::complex<double> > chi2::funcs(double m,std::vector<double> &par);

#ifdef ADOL_ON
template std::vector<std::complex<adouble> > chi2::amps(double m, std::vector<std::complex<adouble> > &cpl,std::vector<adouble> &par);
template std::vector<std::complex<adouble> > chi2::funcs(double m,std::vector<adouble> &par);
#endif//ADOL_ON

std::vector<double> chi2::phase_space(double m){ // gives a vector with phase space factors for each wave at m3pi = m
	double global_ps = phaseSpace(m,_globalPs,0,0.);
	std::vector<double> ps = std::vector<double>(_nWaves);
	for(int wave=0;wave<_nWaves;wave++){
		ps[wave]=global_ps*phaseSpace(m,_wavePs[wave],0,0.);
	};
	return ps;
};


std::string className(){
	return "chi2";
};

void chi2::add_wave(){ // Adds a wave, increases the internal definitions
	_nWaves+=1;
	unsigned int nBor = _borders_waves.size();
	if (nBor==0){
		_borders_waves.push_back(0);
	}else{
		_borders_waves.push_back(_borders_waves[nBor-1]);
	};
	_L.push_back(DEFAULT_L);
	_upperLims.push_back(UPPER_MASS_LIMIT);
	_lowerLims.push_back(LOWER_MASS_LIMIT);
	_wavePs.push_back(0);
	_waveNames.push_back("unnamed_wave");
};
void chi2::add_func(int i){ // Adds a function -  Sets internal definitions accordingly
	_nFuncs+=1;
	_funcs.push_back(i);// Add new
	_L_func.push_back(DEFAULT_L);
	int nPar = getNpars(i);
	if (_maxNpars < nPar){
		_maxNpars = nPar;
	};
	int nParNon=getNparsNonConst(i);
	unsigned int nBor = _borders_par.size();
	if (nBor == 0){
		_borders_par.push_back(nParNon);
	}else{
		int nUp = _borders_par[nBor-1];
//		std::cout<<nUp<<"::::"<<nParNon<<std::endl;
		_borders_par.push_back(nUp+nParNon);
	};
	nBor = _borders_const.size();
	if (nBor == 0){
		_borders_const.push_back(nPar-nParNon);
	}else{
		int nUp = _borders_const[nBor-1];
		_borders_const.push_back(nUp+nPar - nParNon);
	};
	_funcNames.push_back("unnamed_function");
	_funcUpperLims.push_back(UPPER_MASS_LIMIT);
	_funcLowerLims.push_back(LOWER_MASS_LIMIT);
	for (int pn=0; pn<nParNon; pn++){
		_parNames.push_back("unnamed_parameter");
	};
	for (int pn=0;pn<nPar-nParNon;pn++){
		_constNames.push_back("unnamed_constant");
		_const.push_back(0.);
	};
};

void chi2::add_func_to_wave(int wave, int func){ // Sets func to be employed by wave. Set internal definitions accordingly
	int border = _borders_waves[wave];
	std::vector<int> new_relations; // Make new relations
	unsigned int nRel = _funcs_to_waves.size(); // Number of relations up to now
	for (unsigned int i=0; i<nRel;i++){
		new_relations.push_back(_funcs_to_waves[i]); // Up to the end of the block for wave, reproduce old relations
		if (i==border-1){
			new_relations.push_back(func); // Add the new relations at the right point
		};
	};
	if (new_relations.size() ==0){ // If non were set yet, start with one
		new_relations.push_back(func);
	};
	_funcs_to_waves=new_relations; // Set new relations
	for (int i=wave;i<_nWaves;i++){ // INcrease block size for 'wave'
		_borders_waves[i]+=1;
	};
	updateFuncLims(); // Update internal definitions
	updateFuncSpin();
	updateNftw();
};

std::vector<int> chi2::get_wave_functions(int wave){ // Gets all functions used by wave
	std::vector<int> wafe;
	if(_borders_waves.size()==0){
		return wafe;
	};
	int upper = _borders_waves[wave];
	int lower=0;
	if (wave>0){
		lower = _borders_waves[wave-1];
	};
	for (int nnn = lower;nnn<upper;nnn++){
		int coil = _funcs_to_waves[nnn];
		bool already = false;
		for (int aaa=0;aaa<wafe.size();aaa++){
			if (wafe[aaa] == coil){
				already = true;
				break;
			};
		};
		if (not already){
			wafe.push_back(coil);	
		};
	};	
	return wafe;
};
std::vector<int> chi2::get_wave_pars(int wave){ // //  // Gets the parameters numbers for wave
	std::vector<int> pars;
	std::vector<int> funcs = get_wave_functions(wave);
	unsigned int nFu = funcs.size();
	for (unsigned int i=0;i<nFu;i++){
		std::vector<int> funcPars = get_function_pars(funcs[i]);
		pars.insert(pars.end(),funcPars.begin(),funcPars.end());
	};
	return pars;
};
std::vector<int> chi2::get_wave_const(int wave){ //  // Gets the constants numbers for wave
	std::vector<int> consts;
	std::vector<int> funcs = get_wave_functions(wave);
	unsigned int nFu = funcs.size();
	for (unsigned int i=0;i<nFu;i++){
		std::vector<int> funcConst = get_function_const(funcs[i]);
		consts.insert(consts.end(),funcConst.begin(),funcConst.end());
	};
	return consts;
};
std::vector<int> chi2::get_function_pars(int func){ // Gets the paramters numbers for func
	std::vector<int> pars;
	if (_borders_par.size()==0){
		return pars;
	};
	int upper = _borders_par[func];
	int lower = 0;
	if (func >0){
		lower = _borders_par[func-1];
	};
	for (int nnn = lower;nnn<upper;nnn++){
		pars.push_back(nnn);
	};
	return pars;
};
std::vector<int> chi2::get_nPars(){ // gets the numbers of parameters for each func
	std::vector<int> nPars;
	nPars.push_back(_borders_par[0]);
	unsigned int i_max = _borders_par.size();
	for (unsigned int i=1;i<i_max;i++){	
		int diff = _borders_par[i]-_borders_par[i-1];
		nPars.push_back(diff);
	};
	return nPars;
};
std::vector<int> chi2::get_nConst(){ // gets the numbers of constants for each func
	std::vector<int> nConst;
	nConst.push_back(_borders_const[0]);
	int i_max = _borders_const.size();
	for (unsigned int i=1;i<i_max;i++){
		int diff = _borders_const[i]-_borders_const[i-1];
		nConst.push_back(diff);
	};
	return nConst;
};
std::vector<int> chi2::get_function_const(int func){ // Get the constant numbers for func
	std::vector<int> consts;
	if (_borders_const.size()==0){
		return consts;
	};
	int upper = _borders_const[func];
	int lower = 0;
	if (func >0){
		lower = _borders_const[func-1];
	};
	for (int nnn=lower;nnn<upper;nnn++){
		consts.push_back(nnn);
	};
	return consts;
};
std::vector<int> chi2::get_function_waves(int func){ // Get the waves that aemploy the function number 'func'
	std::vector<int> waves;
	for (int wave=0;wave<_nWaves;wave++){
		std::vector<int> funcs = get_wave_functions(wave);
		int i_max = funcs.size();
		for (unsigned int i=0;i<i_max;i++){
			if (funcs[i] == func){
				waves.push_back(wave);
			};
		};
	};
	return waves;
};

void chi2::setWaveName(int i,std::string name){ // Simple setter
	_waveNames[i] = name;
};
void chi2::setFunctionName(int i, std::string name){ // Simple setter
	_funcNames[i] = name;
};
void chi2::setParameterName(int i, std::string name){ // Simple setter
	_parNames[i] = name;
};
void chi2::setConstantName(int i, std::string name){ // Simple setter
	_constNames[i] = name;
};
std::string chi2::getWaveName(int i){ // Simple getter
	return _waveNames[i];
};
std::string chi2::getFunctionName(int i){ // Simple getter
	return _funcNames[i];
};
std::string chi2::getParameterName(int i){ // Simple getter
	return _parNames[i];
};
std::string chi2::getConstantName(int i){ // Simple getter
	return _constNames[i];
};
void chi2::setWaveLimits(int i, double lower, double upper){ // Simple setter
	_upperLims[i] = upper;
	_lowerLims[i] = lower;
	updateFuncLims();
};
void chi2::setWaveSpin(int i, int L){// Simple setter
	_L[i]=L;
	updateFuncSpin();
};
void chi2::setGlobalPhaseSpace(int i){// Simple setter
	_globalPs=i;
};
void chi2::setWavePhaseSpace(int i, int ps){// Simple setter
	_wavePs[i] = ps;
};
void chi2::setConst(int i,double con){// Simple setter
	_const[i] = con;
};

void chi2::updateFuncLims(){ // Updates the mass limits for the funtions from the mass limits for the waves 
	for (int func=0;func<_nFuncs;func++){
			double upper =  LOWER_MASS_LIMIT;// This is exchanged (upper <-> LOWER), since
			double lower =  UPPER_MASS_LIMIT;// the smallest mass region is wanted
			std::vector<int> waves = get_function_waves(func);
			int i_max = waves.size();
			for (unsigned int i=0;i<i_max;i++){
				int wave = waves[i];
				double waveUpper = _upperLims[wave];
				double waveLower = _lowerLims[wave];
				if (waveUpper > upper){
					upper = waveUpper;
				};
				if (waveLower < lower){
					lower = waveLower;
				};
			};
		_funcUpperLims[func] = upper;
		_funcLowerLims[func] = lower;
	};
};

void chi2::updateFuncSpin(){ // Updates the function spin according to their corresponding wave spin (if there is a mistake (different spins for the same function), nothing will happen)
	for (int func=0;func<_nFuncs;func++){
		std::vector<int> waves = get_function_waves(func);
		int i_max = waves.size();
		if (i_max >0){
			int L = _L[waves[0]];
//			for (int i=1;i<i_max;i++){
//				if (_L[waves[i]] != L){
//					std::cout << "Warning: Spins for the same function differ for different waves: _L["<< waves[0] <<"] = "<<_L[waves[0]] <<" != _L[" << waves[i] << "] = "<<_L[waves[i]]<<std::endl;
//				};
//			};
			_L_func[func] = L;
//		}else{
//			std::cout <<"No wave for function "<<func<< " declared. No spin set." << std::endl;
		};
	};
};

bool chi2::checkConsistency(){ // Checks for internal consistency
	int nErr =0;
	if(_nWaves != _waveNames.size()){
		nErr+=1;
		std::cout<< "Inconsistency found: _nWaves does not match _waveNames.size()"<<std::endl;
	};
	if(_nFuncs != _funcNames.size()){
		nErr+=1;
		std::cout << "Inconsistency found: _nFuncs does not match _funcNames.size()"<<std::endl;
	};
	if (_parNames.size() >0){
		int maxPar = _borders_par.back();
		if (_parNames.size() != maxPar){
			nErr+=1;
			std::cout << "Inconsistency found: _parNames.size() does not match _borders_par.back()"<<std::endl;
		};	
	};
	if (_constNames.size() !=0){
		int maxConst = _borders_const.back();
		if(_constNames.size() != maxConst){
			nErr+=1;
			std::cout << "Inconsistency found: _constNames.size() does not match _borders_const.back()"<<std::endl;
		};

	};
	if( _nFuncs != _funcs.size()){
		nErr+=1;
		std::cout << "Inconsistecy found: _nFuncs = "<<_nFuncs<<" does not match _funcs.size() = "<<_funcs.size()<<std::endl;
	};
	for (unsigned int i=1;i<_borders_waves.size();i++){
		if (_borders_waves[i] < _borders_waves[i-1]){
			nErr+=1;
			std::cout << "Inconsistency found: _borders_waves[] is not sorted"<<std::endl;
		};
	};
	for (unsigned int i=1;i<_borders_par.size();i++){
		if (_borders_par[i] < _borders_par[i-1]){
			nErr+=1;
			std::cout << "Inconsistency found: _borders_par[] is not sorted"<<std::endl;
		};
	};
	for (unsigned int i=1;i<_borders_const.size();i++){
		if (_borders_const[i] < _borders_const[i-1]){
			nErr+=1;
			std::cout << "Inconsistency found: _borders_const[] is not sorted"<<std::endl;
		};
	};
	for (unsigned int i=0;i<_funcs_to_waves.size();i++){
		if (_funcs_to_waves[i] > _nFuncs-1){
			nErr+=1;	
			std::cout << "Inconsistency found: Function with number "<<_funcs_to_waves[i]<< " does not exist"<<std::endl;
		};
	};
	if (_wavePs.size() != _nWaves){
		nErr+=1;
		std::cout << "Inconsistency found: _wavePs.size() does not match _nWaves" << std::endl;
	};
	if(_nWaves != _L.size()){
		nErr+=1;
		std::cout<< "Inconsistency found: _nWaves does not match _L.size()"<<std::endl;
	};
	if(_nWaves != _upperLims.size()){
		nErr+=1;
		std::cout<< "Inconsistency found: _nWaves does not match _upperLims.size()"<<std::endl;
	};
	if(_nWaves != _lowerLims.size()){
		nErr+=1;
		std::cout<< "Inconsistency found: _nWaves does not match _lowerLims.size()"<<std::endl;
	};	if (nErr>0){
		return false;
	};
	return true;
};

void chi2::printStatus(){ // Prints the internal status
	std::cout << "_nWaves: " <<_nWaves<<std::endl;
	std::cout << std::endl;
	std::cout << "_nFuncs: " << _nFuncs<<std::endl;
	std::cout << std::endl;
	std::cout<<"_maxNpars: "<<_maxNpars<<std::endl;
	std::cout<<std::endl;
	std::cout<<"_nPar: "<<_nPar<<std::endl;
	std::cout<<std::endl;
	std::cout<<"_nFtw (number of funtions to waves): "<<_nFtw<<std::endl;
	std::cout<<std::endl;
	std::cout << "_funcs"<<std::endl;
	print_vector(_funcs);
	std::cout << std::endl;
	std::cout << "_borders_waves" << std::endl;
	print_vector(_borders_waves);
	std::cout << std::endl;
	std::cout << "_funcs_to_waves" << std::endl;
	print_vector(_funcs_to_waves);
	std::cout << std::endl;
	std::cout << "_borders_par" << std::endl;
	print_vector(_borders_par);
	std::cout << std::endl;
	std::cout << "_borders_const" << std::endl;
	print_vector(_borders_const);
	std::cout << std::endl;
	std::cout << "_L" << std::endl;
	print_vector(_L);
	std::cout << std::endl;
	std::cout<<"_L_func"<<std::endl;
	print_vector(_L_func);
	std::cout<<std::endl;
	std::cout << "_upperLims" << std::endl;
	print_vector(_upperLims);
	std::cout << std::endl;
	std::cout << "_lowerLims" << std::endl;
	print_vector(_lowerLims);
	std::cout << std::endl;
	std::cout << "_funcUpperLims" << std::endl;
	print_vector(_funcUpperLims);
	std::cout << std::endl;
	std::cout << "_funcLowerLims" << std::endl;
	print_vector(_funcLowerLims);
	std::cout << std::endl;
	std::cout << "_const" << std::endl;
	print_vector(_const);
	std::cout << std::endl;	
	std::cout<<"_interface"<<std::endl;
	print_vector(_interface);
	std::cout<<std::endl;
	std::cout<<"_globalPs: "<<_globalPs<<std::endl;
	std::cout << std::endl;
	std::cout << "_wavePs" << std::endl;
	print_vector(_wavePs);
	std::cout << std::endl;
	std::cout << "_waveNames" << std::endl;
	print_vector(_waveNames);
	std::cout << std::endl;
	std::cout << "_funcNames" << std::endl;
	print_vector(_funcNames);
	std::cout << std::endl;
	std::cout << "_parNames" << std::endl;
	print_vector(_parNames);
	std::cout << std::endl;
	std::cout << "_constNames" << std::endl;
	print_vector(_constNames);
	std::cout << std::endl;
	std::cout<<"_write_out: "<<_write_out<<std::endl;
	std::cout<<std::endl;
};

void chi2::printParameters(){ // Prints the paramters in a 'nice' way
	int func_count = 0;
	for (int i=0;i<_nWaves;i++){
		std::cout<<"────────────────────────────────────────────────"<<std::endl;
		std::cout<< _waveNames[i]<<std::endl;
		std::cout<<"┬───────────────────────────────────────────────"<<std::endl;
		std::vector<int> funcs = get_wave_functions(i);
		int j_max = funcs.size();
		for (unsigned int j=0;j<j_max;j++){
			int func = funcs[j];
			if (j==funcs.size()-1){
				std::cout << "└─";
			}else{
				std::cout << "├─";
			};
			std::cout << _funcNames[func]<<" #"<< func_count  << std::endl;
			func_count++;
			std::vector<int> pars = get_function_pars(func);
			std::vector<int> cons = get_function_const(func);
			int k_max = pars.size();
			for(unsigned int k=0;k<k_max;k++){
				int par = pars[k];
				if (j==funcs.size()-1){
					std::cout<<"  ";
				}else{
					std::cout <<"│ ";
				};
				if (k==pars.size()-1 and cons.size()==0){
					std::cout <<"└─";
				}else{
					std::cout <<"├─";
				};
				std::cout<< _parNames[par]<<std::endl;
			};
			k_max = cons.size();
			for (unsigned int k=0;k<k_max;k++){
				int con = cons[k];
				if (j==funcs.size()-1){
					std::cout<<"  ";
				}else{
					std::cout <<"│ ";
				};
				if (k==cons.size()-1){
					std::cout <<"└─";
				}else{
					std::cout <<"├─";
				};
				std::cout<< _constNames[con]<<std::endl;
			};

		};
	std::cout<<std::endl;
	};

};

void chi2::setInterface(){ // Sets up the interface for the 'old' paramter-vector format
	_interface = std::vector<int>();
	for (int wave = 0;wave<_nWaves;wave++){
		int wave_up = _borders_waves[wave];
		int wave_low =0;
		if (wave >0){
			wave_low = _borders_waves[wave-1];	
		};
		for(int par =wave_low; par < wave_up;par++){
			_interface.push_back(1); // real part of cpl
			_interface.push_back(2); // imag part of cpl
			int func = _funcs_to_waves[par];
			int up_func = _borders_par[func];
			int up_con = _borders_const[func];
			int low_func = 0;
			int low_con=0;
			if(func>0){
				low_func = _borders_par[func-1];
				low_con = _borders_const[func-1];
			};

			for(int i=low_func; i<up_func; i++){
				_interface.push_back(3); // parameter
			};
			for(int i=low_con;i<up_con;i++){
				_interface.push_back(4); // constant
			};
		};
	};
	_nPar=0;
	for(unsigned int i=0;i<_interface.size();i++){
		if(_interface[i] ==3){
			_nPar+=1;
		};
		if(_interface[i] ==1){
		};
	};
};
std::vector<std::complex<double> > chi2::cpls(std::vector<double> &par){ // gets the couplings from a paramters-vector in the 'old' format
	std::vector<std::complex<double> > cpl;
	for (unsigned int i =0;i<par.size();i++){
		if (_interface[i]==1){
			cpl.push_back(std::complex<double>(par[i],par[i+1]));
		};
	};
	return cpl;
};
std::vector<double> chi2::pars(std::vector<double> &par){ // gets the shape paramters from a paramters-vector in the 'old' format
	std::vector<double> pars;
	for(unsigned int i=0;i<par.size();i++){
		if (_interface[i] == 3){
			pars.push_back(par[i]);
		};
	};
	return pars;
};

void chi2::updateNftw(){ // Updates the number of function-to-wave couplings
	_nFtw = _funcs_to_waves.size();
};

void chi2::open_output(std::string name){
	_outStream = new ofstream;
	_outStream->open(name.c_str());
	_write_out = true;
};
void chi2::close_output(){
	_outStream->close();
	_write_out = false;
};

int chi2::getNftw(){
	return _nFtw;
};


