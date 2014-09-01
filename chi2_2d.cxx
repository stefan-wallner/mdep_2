#include<iostream>
#include"chi2_2d.h"
#include"chi2.h"
#include"../chi_squared/breitWigners.h"
#include"../chi_squared/phaseSpace.h"

double UPPER_MASS_LIMIT = 2.5;
double LOWER_MASS_LIMIT = 0.5;
int DEFAULT_L =0;


template<typename T>
void print_vector(std::vector<T> in){
	if (in.size()==0){
		std::cout << "Empty vector"<<std::endl;
		return;
	};
	std::cout << "[" << in[0] ;
	for (unsigned int i=1;i<in.size();i++){
		std::cout << ", " << in[i];
	};
	std::cout<<"]"<<std::endl;
};


chi2_2d::chi2_2d():chi2(),_nPoints(0),_nIso(0),_maxNparsIso(0),_maxBinIso(0),_has_isobars(false){};

template<typename xdouble>
std::vector<std::complex<xdouble> > chi2_2d::amps(double m,std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par, std::vector<std::vector<std::complex<xdouble> > > &funcEvals2pi){ // the amplitude for each wave and mass bin
	std::vector<std::complex<xdouble> > funcEval = funcs(m,par); // Evalulated BW-functions
	std::vector<std::complex<xdouble> > ampl = std::vector<std::complex<xdouble> >(_nPoints,std::complex<xdouble>(0.,0.)); // Vector of final amplitudes
	std::vector<double> ps = phase_space(m);
	std::complex<xdouble> amp; // Actual wave amplitude
	int upBor=0; // Upper limit for function number
	int loBor=0; // Lower limit for function number
	int amplcount=0;
	for (int wave =0;wave<_nWaves;wave++){
		int n_iso_bin = _wave_binning_pts[wave];
		loBor = upBor;
		upBor = _borders_waves[wave];	
		if (m >= _lowerLims[wave] and m<_upperLims[wave]){ 
			if(-1==n_iso_bin){
				std::complex<xdouble> amp(0.,0.);	
				for (int nFunc = loBor; nFunc<upBor; nFunc++){
					int func = _funcs_to_waves[nFunc]; // Number of function contributing to the actual wave
					amp+=std::complex<xdouble>(ps[wave],0.)*cpl[nFunc]*funcEval[func];
				};
				ampl[amplcount]=amp;
				amplcount++;
			}else{
				for(int bin2=0;bin2<n_iso_bin;bin2++){
					std::complex<xdouble> amp(0.,0.);
					for (int nFunc = loBor; nFunc<upBor;nFunc++){
						int func =  _funcs_to_waves[nFunc];
						int iso  =  _iso_to_waves[nFunc];
						amp+=std::complex<xdouble>(ps[wave],0.)*cpl[nFunc]*funcEval[func]*funcEvals2pi[iso][bin2];
					};
					ampl[amplcount]=amp;
					amplcount++;
				};
			};
		}else{
				amplcount+=abs(n_iso_bin);
		};
	};
	return ampl;
};
template std::vector<std::complex<double> > chi2_2d::amps(double m, std::vector<std::complex<double> > &cpl,std::vector<double> &par, std::vector<std::vector<std::complex<double> > > &funcEvals2pi);

template<typename xdouble>
std::vector<std::vector<std::complex<xdouble> > > chi2_2d::iso_funcs(std::vector<xdouble> &par){ // Evaluates the isobar paramterizations at ALL masses, so they do not have to be recalculated each time
	if (_has_isobars){
		std::vector<std::vector<std::complex<xdouble> > > f = std::vector<std::vector<std::complex<xdouble> > >(_nIso,std::vector<std::complex<xdouble> >(_maxBinIso,std::complex<xdouble>(0.,0.)));
		int upPar=0;
		int loPar=0;
		int upConst=0;
		int loConst=0;
		std::vector<xdouble> act_par = std::vector<xdouble>(_maxNparsIso);
		for (int iso=0;iso<_nIso;iso++){
			loPar = upPar;
			upPar = _iso_borders_par[iso];
			loConst = upConst;
			upConst = _iso_borders_const[iso];
			int nBins = _iso_binning_pts[iso];
			int nBinning = _iso_n_binning[iso];
			int pos=0;
			int nFunc = _isos[iso];
			for (int i= loPar;i<upPar;i++){
				act_par[pos]=par[i];
				pos++;
			};
			for (int i=loConst;i<upConst;i++){
				act_par[pos]=_iso_const[i];
				pos++;
			};
			for(int bin=0;bin<nBins;bin++){
				double m = (_iso_binnings[nBinning][bin]+_iso_binnings[nBinning][bin+1])/2;
				f[iso][bin]=bw(m,act_par,nFunc,_L_iso[iso]);
			};
		};
		return f;
	};
	return std::vector<std::vector<std::complex<xdouble> > >();
};
template std::vector<std::vector<std::complex<double> > > chi2_2d::iso_funcs(std::vector<double> &par);

#ifdef ADOL_ON // Enables autodiff
template std::vector<std::complex<adouble> > chi2_2d::amps(double m, std::vector<std::complex<adouble> > &cpl,std::vector<adouble> &par, std::vector<std::vector<std::complex<adouble> > > &funcEvals2pi);
template std::vector<std::vector<std::complex<adouble> > > chi2_2d::iso_funcs(std::vector<adouble> &par);
#endif//ADOL_ON

std::string className(){
	return "chi2_2d";
};

void chi2_2d::add_wave(){ // Adds a wave
	chi2::add_wave();
	_L_iso.push_back(DEFAULT_L);
	_wave_n_binning.push_back(-1);
	_wave_binning_pts.push_back(-1);
};

void chi2_2d::add_func_to_wave(int wave, int func){ // Adds a function without an isobar to a wave
	int border = _borders_waves[wave];
	std::vector<int> new_relations;
	std::vector<int> new_relations_iso;
	unsigned int nRel = _funcs_to_waves.size();
	for (unsigned int i=0; i<nRel;i++){
		new_relations.push_back(_funcs_to_waves[i]);
		new_relations_iso.push_back(_iso_to_waves[i]);
		if (i==border-1){
			new_relations.push_back(func);
			new_relations_iso.push_back(-1);
		};
	};
	if (new_relations.size() ==0){
		new_relations.push_back(func);
		new_relations_iso.push_back(-1);
	};
	_funcs_to_waves=new_relations;
	_iso_to_waves=new_relations_iso;
	for (int i=wave;i<_nWaves;i++){
		_borders_waves[i]+=1;
	};
	updateFuncLims();
	updateFuncSpin();
	updateIsobar();
	updateNftw();
};
void chi2_2d::add_funcs_to_wave(int wave, int func, int func2){ // Adds a function-isobar pair to a wave
	int border = _borders_waves[wave];
	std::vector<int> new_relations;
	std::vector<int> new_relations_iso;
	unsigned int nRel = _funcs_to_waves.size();
	for (unsigned int i=0; i<nRel;i++){
		new_relations.push_back(_funcs_to_waves[i]);
		new_relations_iso.push_back(_iso_to_waves[i]);
		if (i==border-1){
			new_relations.push_back(func);
			new_relations_iso.push_back(func2);
			_has_isobars=true;
		};
	};
	if (new_relations.size() ==0){
		new_relations.push_back(func);
		new_relations_iso.push_back(func2);
		_has_isobars=true; // Set falg to indicate isobar parametrization
	};
	_funcs_to_waves=new_relations;
	_iso_to_waves=new_relations_iso;
	for (int i=wave;i<_nWaves;i++){
		_borders_waves[i]+=1;
	};
	updateFuncLims();
	updateFuncSpin();
	updateIsobar();
	updateNftw();
};

void chi2_2d::printStatus(){ // Prints the internal status
	chi2::printStatus();
	std::cout<<std::endl<<"_nPoints: "<<_nPoints<<std::endl;
	std::cout<<std::endl<<"_has_isobars: "<<_has_isobars<<std::endl;
	std::cout<<std::endl<<"_nIso: "<<_nIso<<std::endl;
	std::cout<<std::endl<<"_maxNparsIso: "<<_maxNparsIso<<"; _maxBinIso: "<<_maxBinIso<<std::endl;
	std::cout<<std::endl<<"_isos"<<std::endl;
	print_vector(_isos);
	std::cout<<std::endl<<"_iso_to_waves"<<std::endl;
	print_vector(_iso_to_waves);
	std::cout<<std::endl<<"_iso_borders_par"<<std::endl;
	print_vector(_iso_borders_par);
	std::cout<<std::endl<<"_iso_borders_const"<<std::endl;
	print_vector(_iso_borders_const);
	std::cout<<std::endl<<"_L_iso"<<std::endl;
	print_vector(_L_iso);
	std::cout<<std::endl<<"_L_iso_func"<<std::endl;
	print_vector(_L_iso_func);
	std::cout<<std::endl<<"_iso_const"<<std::endl;
	print_vector(_iso_const);
	std::cout<<std::endl<<"_iso_binning_pts"<<std::endl;
	print_vector(_iso_binning_pts);
	std::cout<<std::endl<<"_wave_binning_pts"<<std::endl;
	print_vector(_wave_binning_pts);
	std::cout<<std::endl<<"_point_to_wave"<<std::endl;
	print_vector(_point_to_wave);
	std::cout<<std::endl<<"_iso_n_binning"<<std::endl;
	print_vector(_iso_n_binning);
	std::cout<<std::endl<<"_wave_n_binning"<<std::endl;
	print_vector(_wave_n_binning);
	std::cout<<std::endl<<"_iso_binnings"<<std::endl;
	for (unsigned int i=0;i<_iso_binnings.size();i++){
		print_vector(_iso_binnings[i]);
	};
	std::cout<<std::endl<<"_point_borders_wave"<<std::endl;
	print_vector(_point_borders_wave);
	std::cout<<std::endl<<"_iso_funcNames"<<std::endl;
	print_vector(_iso_funcNames);
	std::cout<<std::endl<<"_iso_parNames"<<std::endl;
	print_vector(_iso_parNames);
	std::cout<<std::endl<<"_iso_constNames"<<std::endl;
	print_vector(_iso_constNames);
};

bool chi2_2d::checkConsistency(){ // Checks consistency
	int nErr =0;
	if (not chi2::checkConsistency()){
		nErr++;
	};
	if (_iso_to_waves.size()>0){
		if(_iso_to_waves[0] != -1){
			nErr++;
			std::cout<<"Inconsistency found: First wave is de-isobarred"<<std::endl;
		};
	};
	for (int iso=0;iso<_nIso;iso++){
		std::vector<int> waves = get_isobar_waves(iso);
		if(waves.size()>0){
			int binn = _wave_n_binning[waves[0]];
			for (int wave=1;wave<waves.size();wave++){
				if (_wave_n_binning[waves[wave]] != binn){
					nErr++;
					std::cout<<"Inconsistency found: Different binnings for the same isobar"<<std::endl;
				};
			};
		};
	};
	if (nErr>0){
		return false;
	};
	return true;
};

////////////// 2d specials ////////////////////
int chi2_2d::getNpoints(){ // Getter
	return _nPoints;
};

void chi2_2d::updateNpoints(){ // Updates the number of employed points
	int nnn = 0;
	std::vector<int> point_to_wave;
	_point_borders_wave=std::vector<int>(0);
	for (int wave=0;wave<_nWaves;wave++){
		int plus = abs(_wave_binning_pts[wave]);
		nnn+=plus;
		_point_borders_wave.push_back(nnn);
		for (int i=0;i<plus;i++){
			point_to_wave.push_back(wave);
		};
//		std::cout << "le incrise: "<<abs(_wave_binning_pts[wave])<<std::endl;
	};
	_nPoints = nnn;
	_point_to_wave = point_to_wave;
};

void chi2_2d::add_isobar_binning(std::vector<double> binning){ // Adds a new isobar binning
	_iso_binnings.push_back(binning);
	int nnew = binning.size()-1;
	if (_maxBinIso<nnew){
		_maxBinIso=nnew;
	};
};
void chi2_2d::setWaveIsobarSpin(int wave, int L){ // Sets the isobar spin for a wave
	_L_iso[wave]=L;
	updateIsobar();
};
void chi2_2d::add_iso(int i){ // Add a new isobar, sets the internal definitions accordingly
	_nIso+=1;
	_isos.push_back(i);
	int nPar = getNpars(i);
	if (_maxNparsIso < nPar){
		_maxNparsIso = nPar;
	};
	int nParNon=getNparsNonConst(i);
	unsigned int nBor = _iso_borders_par.size();
	if (nBor == 0){
		_iso_borders_par.push_back(nParNon);
	}else{
		int nUp = _iso_borders_par[nBor-1];
//		std::cout<<nUp<<"::::"<<nParNon<<std::endl;
		_iso_borders_par.push_back(nUp+nParNon);
	};
	nBor = _iso_borders_const.size();
	if (nBor == 0){
		_iso_borders_const.push_back(nPar-nParNon);
	}else{
		int nUp = _iso_borders_const[nBor-1];
		_iso_borders_const.push_back(nUp+nPar - nParNon);
	};
	_iso_funcNames.push_back("unnamed_function");
	for (int pn=0; pn<nParNon; pn++){
		_iso_parNames.push_back("unnamed_parameter");
	};
	for (int pn=0;pn<nPar-nParNon;pn++){
		_iso_constNames.push_back("unnamed_constant");
		_iso_const.push_back(0.);
	};
	_L_iso_func.push_back(DEFAULT_L);
	_iso_binning_pts.push_back(-1);
	_iso_n_binning.push_back(-1);
};

std::vector<int> chi2_2d::get_nParsIso(){ // gets the numbers of paramters for each isobar
	std::vector<int> nPars;
	nPars.push_back(_iso_borders_par[0]);
	unsigned int i_max = _iso_borders_par.size();
	for (unsigned int i=1;i<i_max;i++){	
		int diff = _iso_borders_par[i]-_iso_borders_par[i-1];
		nPars.push_back(diff);
	};
	return nPars;
};
std::vector<int> chi2_2d::get_nConstIso(){ // gets the numbers of constants for each isobar
	std::vector<int> nConst;
	nConst.push_back(_iso_borders_const[0]);
	int i_max = _iso_borders_const.size();
	for (unsigned int i=1;i<i_max;i++){
		int diff = _iso_borders_const[i]-_iso_borders_const[i-1];
		nConst.push_back(diff);
	};
	return nConst;
};

std::vector<int> chi2_2d::get_wave_isobars(int wave){ // Gets all isobars in wave
	int upper = _borders_waves[wave];
	int lower=0;
	if (wave>0){
		lower = _borders_waves[wave-1];
	};
	std::vector<int> isos;
	for (int i=lower;i<upper;i++){
		bool already = false;
		for (int aaa=0;aaa<isos.size();aaa++){
			if (isos[aaa]==_iso_to_waves[i]){
				already = false;
				break;
			};
		};
		if (_iso_to_waves[i] >-1 and not already){
			isos.push_back(_iso_to_waves[i]);
		};
	};
	return isos;
};
std::vector<int> chi2_2d::get_wave_iso_pars(int wave){ // Gets the isobar paramters for wave
	std::vector<int> isos = get_wave_isobars(wave);
	std::vector<int> pars;
	for (int i=0;i<isos.size();i++){
		std::vector<int> isopars = get_isobar_pars(isos[i]);
		pars.insert(pars.end(),isopars.begin(),isopars.end());
	};
	return pars;
};

std::vector<int> chi2_2d::get_wave_iso_const(int wave){ // Gets the isobar constants for wave
	std::vector<int> isos = get_wave_isobars(wave);
	std::vector<int> consts;
	for (int i=0;i<isos.size();i++){
		std::vector<int> isoconst = get_isobar_const(isos[i]);
		consts.insert(consts.end(),isoconst.begin(),isoconst.end());
	};
	return consts;
};

std::vector<int> chi2_2d::get_isobar_pars(int func){ // Gets the numbers of paramters for the 'func' isobar paramterization
	std::vector<int> pars;
	if(0==_iso_borders_par.size()){
		return pars;
	};
	int upper = _iso_borders_par[func];
	int lower = 0;
	if (func>0){
		lower = _iso_borders_par[func-1];
	};
	for (int i=lower;i<upper;i++){
		pars.push_back(i);
	};
	return pars;
};
std::vector<int> chi2_2d::get_isobar_const(int func){ // Gets the numbers of constants for the 'func' isobar paramterization
	std::vector<int> consts;
	if(0==_iso_borders_const.size()){
		return consts;
	};
	int upper = _iso_borders_const[func];
	int lower =0;
	if (func>0){
		lower = _iso_borders_const[func];
	};
	for (int i=lower;i<upper;i++){
		consts.push_back(i);
	};
	return consts;
};
std::vector<int> chi2_2d::get_isobar_waves(int iso){ // Get the waves that employ a certain isobar
	std::vector<int> waves;
	for (int i=0;i<_nWaves;i++){
		std::vector<int> isos = get_wave_isobars(i);
		for (int ff=0; ff<isos.size();ff++){
			if(iso==isos[ff]){
				waves.push_back(i);
				break;
			};
		};		
	};
	return waves;
};
void chi2_2d::updateIsobar(){ // Updates the internal definitions for the isobar paramterizations
	for (int wave =0;wave<_nWaves	;wave++){
		std::vector<int> isos = get_wave_isobars(wave);
		for (int niso=0;niso<isos.size();niso++){
			int iso = isos[niso];
			_L_iso_func[iso]     =_L_iso[wave];
			_iso_binning_pts[iso]=_wave_binning_pts[wave];
			_iso_n_binning[iso]  =_wave_n_binning[wave];
		};
	};
	updateNpoints();
};
void chi2_2d::setIsobarName(int i,std::string name){ // Simple setter
	_iso_funcNames[i]=name;
};
void chi2_2d::setIsoParName(int i,std::string name){ // Simple setter
	 _iso_parNames[i]=name;
};
void chi2_2d::setIsoConstName(int i,std::string name){ // Simple setter
	_iso_constNames[i]=name;
};
std::string chi2_2d::getIsobarName(int i){ // Simple getter
	return _iso_funcNames[i];
};
std::string chi2_2d::getIsoParName(int i){ // Simple getter
	return  _iso_parNames[i];
};
std::string chi2_2d::getIsoConstName(int i){ // Simple getter
	return _iso_constNames[i];
};

void chi2_2d::setWaveIsobarBinning(int wave, int binning){ // Sets the isobar binning for a certain wave
	_wave_n_binning[wave] = binning;
	int npoints = _iso_binnings[binning].size()-1;
	_wave_binning_pts[wave] = npoints;
	updateIsobar();
};

void chi2_2d::set_iso_const(int con, double value){ // Sets constants in isobar paramterizations
	_iso_const[con] = value;
};

void chi2_2d::printParameters(){ // Prints the parameters in a nice way
	int func_count = 0;
	int iso_count =0;
	for (int i=0;i<_nWaves;i++){
		std::cout<<"────────────────────────────────────────────────"<<std::endl;
		std::cout<< _waveNames[i]<<std::endl;
		std::cout<<"┬───────────────────────────────────────────────"<<std::endl;
		std::vector<int> funcs = get_wave_functions(i);
		std::vector<int> isos  = get_wave_isobars(i);
		int j_max = funcs.size();
		int is_max= isos.size();
		for (unsigned int j=0;j<j_max;j++){
			int func = funcs[j];
			if (j==funcs.size()-1 and 0==is_max){
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
				if (j==funcs.size()-1 and 0==is_max){
					std::cout<<"  ";
				}else{
					std::cout <<"│ ";
				};
				if (k==pars.size()-1 and cons.size()==0 ){
					std::cout <<"└─";
				}else{
					std::cout <<"├─";
				};
				std::cout<< _parNames[par]<<std::endl;
			};
			k_max = cons.size();
			for (unsigned int k=0;k<k_max;k++){
				int con = cons[k];
				if (j==funcs.size()-1 and 0==is_max){
					std::cout<<"  ";
				}else{
					std::cout <<"│ ";
				};
				if (k==cons.size()-1 ){
					std::cout <<"└─";
				}else{
					std::cout <<"├─";
				};
				std::cout<< _constNames[con]<<std::endl;
			};

		};
		if (is_max > 0){
			std::cout<<"│ISOBARS:"<<std::endl;
			for (int is=0;is<is_max;is++){
				int iso = isos[is];
				if ( isos.size()-1==is){
					std::cout<<"└─";
				}else{
					std::cout << "├─";
				};
				std::cout << _iso_funcNames[iso]<<" #"<< iso_count  << std::endl;
				iso_count++;
				std::vector<int> pars = get_isobar_pars(iso);
				std::vector<int> cons = get_isobar_const(iso);
				int k_max = pars.size();
				for(unsigned int k=0;k<k_max;k++){
					int par = pars[k];
					if (is==isos.size()-1){
						std::cout<<"  ";
					}else{
						std::cout <<"│ ";
					};
					if (k==pars.size()-1 and cons.size()==0){
						std::cout <<"└─";
					}else{
						std::cout <<"├─";
					};
					std::cout<< _iso_parNames[par]<<std::endl;
				};
				k_max = cons.size();
				for (unsigned int k=0;k<k_max;k++){
					int con = cons[k];
					if (is==isos.size()-1){
						std::cout<<"  ";
					}else{
						std::cout <<"│ ";
					};
					if (k==cons.size()-1){
						std::cout <<"└─";
					}else{
						std::cout <<"├─";
					};
					std::cout<< _iso_constNames[con]<<std::endl;
				};
	
			};
		};
	std::cout<<std::endl;
	};

};


