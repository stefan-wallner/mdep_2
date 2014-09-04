#include<iostream>
#include"waveset.h"
#include"../chi_squared/breitWigners.h"
#include"../chi_squared/phaseSpace.h"

#include"invert33.h"


double UPPER_MASS_LIMIT = 2.5;
double LOWER_MASS_LIMIT = 0.5;
int DEFAULT_L =0;

waveset::waveset():
	_nWaves(0),
	_nFuncs(0),
	_globalPs(0),
	_maxNpars(0),
	_nPar(0),
	_nFtw(0),
	_nPoints(0),
	_nIso(0),
	_maxNparsIso(0),
	_maxBinIso(0),
	_nBranch(0),
	_nTbin(1),
	_nBins(1),
	_write_out(false),
	_has_isobars(false),
	_binning(std::vector<double>(2)){
	setTbinning(std::vector<double>(2,0.));
};
//########################################################################################################################################################
#ifdef USE_YAML
///Constructror from YAML files
waveset::waveset(
							std::string 						card):
	_nWaves(0),
	_nFuncs(0),
	_globalPs(0),
	_maxNpars(0),
	_nPar(0),
	_nFtw(0),
	_nPoints(0),
	_nIso(0),
	_maxNparsIso(0),
	_maxBinIso(0),
	_nBranch(0),
	_nTbin(1),
	_nBins(1),
	_write_out(false),
	_has_isobars(false),
	_binning(std::vector<double>(2)){

	setTbinning(std::vector<double>(2,0.));
	YAML::Node Ycard   = YAML::LoadFile(card);
	std::string parametrizations = Ycard["parametrization_file"].as<std::string>();
	std::string waves = Ycard["wave_file"].as<std::string>();
	YAML::Node Yparam  = YAML::LoadFile(parametrizations);
	YAML::Node Ywaves  = YAML::LoadFile(waves);
	loadGlobalPhaseSpace(Ycard);
	std::map<std::string,int> fMap = loadFunctions(Ycard, Yparam);
	loadWaves(Ycard, Ywaves);
	loadFtw(Ycard, fMap);
	loadBranchings(Ycard);
	loadBinnings(Ycard);
	if(_has_isobars){
		std::string binning_file;
		if (Ycard["isobar_binnings"]){
			binning_file = Ycard["isobar_binnings"].as<std::string>();
		}else{
			std::cerr<<"Error: No isobar binnings given"<<std::endl;
		};
		std::cout<<"open: "<<binning_file<<std::endl;
		YAML::Node Yisob = YAML::LoadFile(binning_file);
		loadIsoBinnings(Ycard,Yisob);
	};
};
#endif//USE_YAML
//########################################################################################################################################################
///Gives the amplitudes for all waves and isobar-mass bins
template<typename xdouble>
std::vector<std::complex<xdouble> > waveset::amps(
							double 							m,
							std::vector<std::complex<xdouble> > 			&cpl,
							std::vector<xdouble> 					&par,
							std::vector<std::vector<std::complex<xdouble> > > 	&funcEvals2pi){

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
template std::vector<std::complex<double> > waveset::amps(double m, std::vector<std::complex<double> > &cpl,std::vector<double> &par, std::vector<std::vector<std::complex<double> > > &funcEvals2pi);
//########################################################################################################################################################
///Returns the function values at m3pi = m and shape parameters par
template<typename xdouble>
std::vector<std::complex<xdouble> > waveset::funcs(
							double 							m,
							std::vector<xdouble> 					&par){

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
template std::vector<std::complex<double> > waveset::funcs(double m,std::vector<double> &par);
//########################################################################################################################################################
///Evaluates the isobar paramterizations at ALL masses, so they do not have to be recalculated each time
template<typename xdouble>
std::vector<std::vector<std::complex<xdouble> > > waveset::iso_funcs(
							std::vector<xdouble> 					&par){

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
template std::vector<std::vector<std::complex<double> > > waveset::iso_funcs(std::vector<double> &par);
//########################################################################################################################################################
///Gives a vector with phase space factors for each wave at m3pi = m
std::vector<double> waveset::phase_space(
								double 							m){

	double global_ps = phaseSpace(m,_globalPs,0,0.);
	std::vector<double> ps = std::vector<double>(_nWaves);
	for(int wave=0;wave<_nWaves;wave++){
		ps[wave]=global_ps*phaseSpace(m,_wavePs[wave],0,0.);
	};
	return ps;
};
//########################################################################################################################################################
///Enable Auto-diff, if needed
#ifdef ADOL_ON
template std::vector<std::complex<adouble> > waveset::amps(double m, std::vector<std::complex<adouble> > &cpl,std::vector<adouble> &par, std::vector<std::vector<std::complex<adouble> > > &funcEvals2pi);
template std::vector<std::complex<adouble> > waveset::funcs(double m,std::vector<adouble> &par);
template std::vector<std::vector<std::complex<adouble> > > waveset::iso_funcs(std::vector<adouble> &par);
#endif//ADOL_ON
//########################################################################################################################################################
///Adds a wave, increases the internal definitions
void waveset::add_wave(){

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
	_L_iso.push_back(DEFAULT_L);
	_wave_n_binning.push_back(-1);
	_wave_binning_pts.push_back(-1);
};
//########################################################################################################################################################
///Adds a function - Sets internal definitions accordingly
void waveset::add_func(
							int 							i, 	// # of function
							bool							is_t_dep){

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
	if (is_t_dep){
		_const_is_t.push_back(_borders_const[_borders_const.size()-1]-1);
	};
};
//########################################################################################################################################################
///Adds a new isobar, sets the internal definitions accordingly
void waveset::add_iso(
							int 							i){ 	// # of function for isobar
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
//########################################################################################################################################################
///Adds a function without an isobar to a wave
void waveset::add_func_to_wave(
							int 							wave,
							int 							func){

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
	if (0==new_relations.size()){
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
	handle_branchings(wave,func);
};
//########################################################################################################################################################
///Adds a function-isobar pair to a wave
void waveset::add_funcs_to_wave(
							int 							wave,
							int 							func,
							int 							func2){ // func2 == isobar

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
	if (0==new_relations.size()){
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
	handle_branchings(wave,func);
};
//########################################################################################################################################################
///Couples two functions: {C1(t), C2(t)} -> {B1*C(t), B2*C(t)}
void waveset::couple_funcs(
							int 						i1,	// # of wave 1
							int 						i2){ 	// # of wave 2

	if (i1 == i2){
		std::cerr<<"Error: Trying to couple two times the same functions"<<std::endl;
		return;
	};
	int branch1 = _coupled[i1];
	int branch2 = _coupled[i2];
	if (branch1 == branch2 and branch1 != -1){
		return;
	};
	if (-1==branch1 and -1==branch2){
		int last_branch = -1;
		int nCpl = _coupled.size();
		for (int i = 0;i<nCpl;i++){
			if (_coupled[i] > last_branch){
				last_branch = _coupled[i];
			};
		};
		last_branch++;
		_coupled[i1] = last_branch;
		_coupled[i2] = last_branch;
	}else if(-1==branch1){
		_coupled[i1] = branch2;
	}else if(-1==branch2){
		_coupled[i2] = branch1;
	}else{
		std::cerr<<"function["<<i1<<"] and function["<<i2<<"] are already in different branchings"<<std::endl;
	};
	update_n_cpls();
	update_n_branch();
};
//########################################################################################################################################################
///Simple setter for the mass limits of the single waves
void waveset::setWaveLimits(
							int 							i, 	// # of wave
							double 							lower,
							double 							upper){

	_upperLims[i] = upper;
	_lowerLims[i] = lower;
	updateFuncLims();
	update_min_max_bin();
};
//########################################################################################################################################################
///Simple setter for the wave spins
void waveset::setWaveSpin(
							int 							i, 	// # of wave
							int 							L){
	_L[i]=L;
	updateFuncSpin();
};
//########################################################################################################################################################
///Sets the isobar spin for a wave
void waveset::setWaveIsobarSpin(
							int 							wave,
							int 							L){

	_L_iso[wave]=L;
	updateIsobar();
};
//########################################################################################################################################################
///Simple setter for the global phase space factor
void waveset::setGlobalPhaseSpace(
							int 							i){	// # of phase space function

	_globalPs=i;
};
//########################################################################################################################################################
///Simple setter for the wave-specific phase space factors
void waveset::setWavePhaseSpace(
							int 							i, 	// # of wave
							int 							ps){	// # of phase space function

	_wavePs[i] = ps;
};
//########################################################################################################################################################
///Sets the isobar binning for a certain wave
void waveset::setWaveIsobarBinning(
							int 							wave,
							int 							binning){

	_wave_n_binning[wave] = binning;
	int npoints = _iso_binnings[binning].size()-1;
	_wave_binning_pts[wave] = npoints;
	updateIsobar();
};
//########################################################################################################################################################
///Simple setter for the constants used in one function
void waveset::setConst(
							int 							i, 	// # of const
							double 							con){

	_const[i] = con;
	for (unsigned int j=0;j<_const_is_t.size();j++){
		if (i == _const_is_t[j]){
			std::cout<<"Warning: Trying to set _const["<<i<<"] which is t' and will be overwritten"<<std::endl;
		};
	};
};
//########################################################################################################################################################
///Sets constants in isobar paramterizations
void waveset::set_iso_const(
							int 							con, 	// # of const
							double 							value){

	_iso_const[con] = value;
};

//########################################################################################################################################################
///Adds a new isobar binning
void waveset::add_isobar_binning(
							std::vector<double> 					binning){

	_iso_binnings.push_back(binning);
	int nnew = binning.size()-1;
	if (_maxBinIso<nnew){
		_maxBinIso=nnew;
	};
};
//########################################################################################################################################################
///Simple setter for wave names
void waveset::setWaveName(
							int 							i, 	// # of wave
							std::string 						name){

	_waveNames[i] = name;
};
//########################################################################################################################################################
///Simple setter for function name
void waveset::setFunctionName(
							int 							i, 	// # of function
							std::string 						name){

	_funcNames[i] = name;
};
//########################################################################################################################################################
///Simple setter for paremeter name
void waveset::setParameterName(
							int 							i, 	// # of parameter
							std::string 						name){
	_parNames[i] = name;
};
//########################################################################################################################################################
///Simple setter for constant name
void waveset::setConstantName(
							int 							i, 	// # of constant
							std::string 						name){

	_constNames[i] = name;
};
//########################################################################################################################################################
///Simple setter for isobar names
void waveset::setIsobarName(
							int 							i,	// # of isobar
							std::string 						name){

	_iso_funcNames[i]=name;
};
//########################################################################################################################################################
///Simple setter for isopar parameter names
void waveset::setIsoParName(
							int 							i,	// # of isobar parameter
							std::string 						name){

	 _iso_parNames[i]=name;
};
//########################################################################################################################################################
///Simple setter for isobar constant names
void waveset::setIsoConstName(
							int 							i,	// # of isobar constant
							std::string 						name){
	_iso_constNames[i]=name;
};
//########################################################################################################################################################
///Sets the binning in m3pi
void waveset::setBinning(
							std::vector<double> 					binning){

	_binning = binning;
	_nBins = binning.size()-1;
	_mMin = binning[0];
	_mMax = binning[_nBins];
	int minBin = 0;
	int maxBin = _nBins-1;
	update_min_max_bin();
};
//########################################################################################################################################################
///Sets the binning in t'
void waveset::setTbinning(
							std::vector<double> 					binning){

	_t_binning = binning;
	_nTbin = binning.size() - 1;
	if (_eval_tbin.size() != _nTbin){
		_eval_tbin = std::vector<bool>(_nTbin,true);
	};
};
//########################################################################################################################################################
///Switches single t' bins on/off
void waveset::setEvalTbin(
							int 							i,	// # of t' bin
							bool 							flag){
	_eval_tbin[i] = flag;
};
//########################################################################################################################################################
#ifdef USE_YAML
///Loads the phase space from a YAML file
bool waveset::loadGlobalPhaseSpace(
							YAML::Node 						&waveset){

	bool sgps = false;
	if (waveset["standard_phase_space"]){
		if (waveset["standard_phase_space"].as<bool>()){
			setGlobalPhaseSpace(10);
			sgps=true;
		};
	};
	int nPS = waveset["global_phase_space"].size();
	for (int i=0;i<nPS;i++){
		setGlobalPhaseSpace(waveset["global_phase_space"][i].as<int>());
		sgps=true;
	};
	if (not sgps){
		std::cout<<"Warning: No global phase space set" << std::endl;
	};
	return sgps;
};
//########################################################################################################################################################
///Loads the functions from a YAML file
std::map<std::string,int> waveset::loadFunctions(
							YAML::Node						&waveset,
							YAML::Node						&param){

	std::map<std::string,int> fMap;
	int fCount = 1; // Starts as one, since map[nonExistingKey] = 0 -> use this to determine, whether a function exists.
	int iCount = 1;
	int pCount = 0;
	int ipCount= 0;
	int cCount = 0;
	int icCount= 0;
	int nWaves = waveset["waves"].size();
	std::vector<double> pVals;
	std::vector<double> pSteps;
	std::vector<bool> pRels;
	for (int i=0;i<nWaves;i++){
		int nFunc = waveset["waves"][i]["parametrizations"].size();
		for (int j=0;j<nFunc;j++){
			std::string fName;
			std::string iName = "none";
			bool isisobar = false;
			if (waveset["waves"][i]["parametrizations"][0].size() != waveset["waves"][i]["parametrizations"][j].size()){
				std::cerr<<"Error: Isobarred and de-isobarred parametrizations in the same wave"<<std::endl;
			};
			if (0==waveset["waves"][i]["parametrizations"][j].size()){
				fName = waveset["waves"][i]["parametrizations"][j].as<std::string>();
//				std::cout<<fName<<" alone"<<std::endl;
			}else if(2==waveset["waves"][i]["parametrizations"][j].size()) {
				isisobar = true;
				fName = waveset["waves"][i]["parametrizations"][j][0].as<std::string>();
				iName = waveset["waves"][i]["parametrizations"][j][1].as<std::string>();
//				std::cout<<iName<<" to "<<fName<<std::endl;
			};
			if(param[fName]){
				if (0==fMap[fName]){
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
					add_func(nParam,isTdep);
					setFunctionName(fCount-1,fName);
					fCount++;
					int nParDef = param[fName]["parameters"].size();
					if (nParDef == nPar){
						for(int par=0;par<nPar;par++){
							std::string pName = param[fName]["parameters"][par]["name"].as<std::string>();
							setParameterName(pCount,pName);
							pCount++;
						};
					}else{
						std::cerr<<"Error: Number of defined parameters does not match required number for "<<fName<<std::endl;
					};
					int nConDef = param[fName]["constants"].size();
					if(nConDef == nCon or (nConDef == nCon-1 and isTdep)){
						for (int con=0;con<nCon;con++){
							std::string cName = param[fName]["constants"][con]["name"].as<std::string>();
							setConstantName(cCount ,cName);
							if (con < nCon -1 or not isTdep){
								double value  = param[fName]["constants"][con]["value"].as<double>();
								setConst(cCount,value);
							};
							cCount++;
						};
					}else{
						std::cerr<<"Error: Number of defined constants does not match required number for "<<fName<<std::endl;
					};
				};
			}else{
				std::cerr << "Error: '"<<fName<<"' not defined in the parametrization-file"<<std::endl;
			};
			if (isisobar){
				if(param[iName]){
					if (0==fMap[iName]){ // Write isobars also in the fMap
						fMap[iName]=iCount;
						int nParam = param[iName]["parametrization"].as<int>();
						int nPar = getNparsNonConst(nParam);
						int nCon = getNpars(nParam)-getNparsNonConst(nParam);
						if (param[iName]["t_dependent"]){
							if (param[iName]["t_dependent"].as<bool>()){
								std::cerr<<"Error: Tryig to set t' dependent isobar parametrization. Not supportet at the moment"<<std::endl;
							};
						};
						add_iso(nParam);
						setIsobarName(iCount-1,iName);
						iCount++;
						int nParDef = param[iName]["parameters"].size();
						if (nParDef == nPar){
							for (int par=0;par<nPar;par++){
								std::string pName = param[iName]["parameters"][par]["name"].as<std::string>();
								setIsoParName(ipCount,pName);
								ipCount++;
							};
						}else{
							std::cerr<<"Error: Number of defined parameters does not match required number for "<<iName<<std::endl;
						};
						int nConDef = param[iName]["constants"].size();
						if(nConDef == nCon){
							for (int con=0;con<nCon;con++){
								std::string cName = param[iName]["constants"][con]["name"].as<std::string>();
								setIsoConstName(icCount,cName);
								double value = param[iName]["constants"][con]["value"].as<double>();
								set_iso_const(icCount,value);
								icCount++;
							};
						}else{
							std::cerr<<"Error: Number of defined constants does not match required number for "<<iName<<std::endl;
						};
					};
				}else{
					std::cerr << "Error: '"<<iName<<"' not defined in the parametrization-file"<<std::endl;
				};
			};
		};
	};
	return fMap;
};
//########################################################################################################################################################
///Loads the waves from a YAML file
bool waveset::loadWaves(
							YAML::Node						&waveset,
							YAML::Node						&defs){

	bool ookk = true;
	int nWaves = waveset["waves"].size();
	double mmin = waveset["waves"][0]["mmin"].as<double>();
	double mmax = waveset["waves"][0]["mmax"].as<double>();
	for (int i=0;i<nWaves;i++){
		add_wave();
		std::string wName = waveset["waves"][i]["name"].as<std::string>();
		setWaveName(i,wName);
		if (not defs[wName]){
			ookk=false;
			std::cerr<<"Error: "<< wName<<" not defined in the wave-definitions file"<<std::endl;
		};
		setWaveSpin(i,defs[wName]["spin"].as<int>());
		if(defs[wName]["isobar_spin"]){
			setWaveIsobarSpin(i,defs[wName]["isobar_spin"].as<int>());
		};
		double mmin_act = waveset["waves"][i]["mmin"].as<double>();
		double mmax_act = waveset["waves"][i]["mmax"].as<double>();
		if (mmin_act<mmin){
			ookk=false;
			std::cerr << "Error: Lower limit of '"<<waveset["waves"][i]["name"]<<"' below anchor mass limit"<<std::endl;
		};
		if (mmax_act>mmax){
			ookk=false;
			std::cerr << "Error: Upper limit of '"<<waveset["waves"][i]["name"]<<"' above anchor mass limit"<<std::endl;
		};
		setWaveLimits(i,mmin_act,mmax_act);
		setWavePhaseSpace(i,defs[wName]["phase_space"].as<int>());
	};
	return ookk;
};
//########################################################################################################################################################
///Puts funtions and waves together // Does not support isobars at the moment
void waveset::loadFtw(
							YAML::Node						&waveset,
							std::map<std::string,int>				&fMap){

	int nWaves = waveset["waves"].size();
	for (int i=0;i<nWaves;i++){
		int nFunc = waveset["waves"][i]["parametrizations"].size();
		for (int func=0;func<nFunc;func++){
			std::string fName;
			std::string iName;
			bool isisobar = false;
			if (0==waveset["waves"][i]["parametrizations"][func].size()){
				fName = waveset["waves"][i]["parametrizations"][func].as<std::string>();
				int nf = fMap[fName]-1;
				add_func_to_wave(i,nf);
//				std::cout<<fName<<" alone"<<std::endl;
			}else if(2==waveset["waves"][i]["parametrizations"][func].size()) {
				isisobar = true;
				fName = waveset["waves"][i]["parametrizations"][func][0].as<std::string>();
				iName = waveset["waves"][i]["parametrizations"][func][1].as<std::string>();
				int nf = fMap[fName]-1;
				int ni = fMap[iName]-1;
				add_funcs_to_wave(i,nf,ni);
//				std::cout<<iName<<" to "<<fName<<std::endl;
			};
		};
	};
};
//########################################################################################################################################################
///Loads branchings from a YAML file
void waveset::loadBranchings(
							YAML::Node						&waveset){

	int nWaves = waveset["waves"].size();
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
					std::string fnam_d;
					if (0==waveset["waves"][wave]["parametrizations"][func].size()){
						fnam_d = waveset["waves"][wave]["parametrizations"][func].as<std::string>();
					}else if(2==waveset["waves"][wave]["parametrizations"][func].size()){
						fnam_d = waveset["waves"][wave]["parametrizations"][func][0].as<std::string>(); // When a coupled function appear with different isobars in one wave, all of them will be "branched"
					};
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
				couple_funcs(wave_numbers[0],wave_numbers[cplto]);
			};
		};
	};
};
//########################################################################################################################################################
///Loads binnings (m3Pi and t') from YAML file
void waveset::loadBinnings(
							YAML::Node						&waveset){

	double M = waveset["mmin"].as<double>();
	double m_max = waveset["mmax"].as<double>();
	double binwidth = waveset["binwidth"].as<double>();
	std::vector<double> binning;
	binning.push_back(M);
	while(M<m_max){
		M+=binwidth;
		binning.push_back(M);
	};
	setBinning(binning);
	std::vector<double> tbinning;
	if (waveset["t_binning"]){
		int nTbin = waveset["t_binning"].size();
		for (int i=0; i<nTbin;i++){
			tbinning.push_back(waveset["t_binning"][i].as<double>());
		};
		setTbinning(tbinning);
	};
	update_min_max_bin();
};
//########################################################################################################################################################
///Load isobar binnings from YAML file
void waveset::loadIsoBinnings(
							YAML::Node 						&waveset,
							YAML::Node						&binnings){

	std::map<std::string,int> binnings_map;
	int nWaves = waveset["waves"].size();
	int nBinn = 1;// Start with 1 for int maps, to determine, if entry exists
	for (int i=0;i<nWaves;i++){ // Create the binnigns map
		std::string wName = waveset["waves"][i]["name"].as<std::string>();
		int nPar = waveset["waves"][i]["parametrizations"].size();
		if (nPar>0){
			int nPer = waveset["waves"][i]["parametrizations"][0].size();
			if (2==nPer){
				if(waveset["waves"][i]["isobar_binning"]){
					std::string binning_name = waveset["waves"][i]["isobar_binning"].as<std::string>();
					if (0==binnings_map[binning_name]){
						binnings_map[binning_name] = nBinn;
						nBinn++;
						std::vector<double> binning;
						int binning_size = binnings[binning_name].size();
						for (int bin;bin<binning_size;bin++){
							binning.push_back(binnings[binning_name][bin].as<double>());
						};
						add_isobar_binning(binning);
					};
					setWaveIsobarBinning(i,binnings_map[binning_name]-1);
				}else{
					std::cerr<<"Error: De-isobarred wave '"<<wName<<"' does not have an isobar binning"<<std::endl;
				};
			};
		};
	};
	


};
#endif//USE_YAML
//########################################################################################################################################################
///Gets total number of points: sum_{waves} n_{isobar_bins} // If no wave is de-isobarred, returns _nWaves
int waveset::getNpoints(){

	return _nPoints;
};
//########################################################################################################################################################
///Returns the number of couplings 'functions to waves'
int waveset::getNftw(){

	return _nFtw;
};
//########################################################################################################################################################
///Gives the number of m3Pi bins
int waveset::getNbins(){

	return _nBins;
};
//########################################################################################################################################################
///Gives the number of t' bins
int waveset::getNtBin(){

	return _nTbin;
};
//########################################################################################################################################################
///Total number of parameters (including non achor couplings)
int waveset::getNtot(){

	return 2*getNcpl() + getNpar() + 2*getNbra() + getNiso();
};
//########################################################################################################################################################
///Number of total couplings
int waveset::getNcpl(){

	return _nBrCpl * _nTbin;
};
//########################################################################################################################################################
///Returns the number of parameters (without branchings or couplings)
int waveset::getNpar(){

	if (0==_borders_par.size()){
		return 0;
	};
	return _borders_par[_borders_par.size() -1];
};
//########################################################################################################################################################
///Number of sets coupled by branchings
int waveset::getNbra(){

	return _nBranch;
};
//########################################################################################################################################################
///Returns the number of parameters for isobars
int waveset::getNiso(){

	int siz = _iso_borders_par.size();
	if (0==siz){
		return 0;
	};
	return _iso_borders_par[siz-1];
};
//########################################################################################################################################################
///Returns the name of waves
std::string waveset::getWaveName(
							int 							i){ 	// # of wave

	return _waveNames[i];
};
//########################################################################################################################################################
/// Gets all functions used by wave
std::vector<int> waveset::get_wave_functions(
							int 							wave){

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
//########################################################################################################################################################
///Gets the parameters numbers for wave
std::vector<int> waveset::get_wave_pars(
							int 							wave){

	std::vector<int> pars;
	std::vector<int> funcs = get_wave_functions(wave);
	unsigned int nFu = funcs.size();
	for (unsigned int i=0;i<nFu;i++){
		std::vector<int> funcPars = get_function_pars(funcs[i]);
		pars.insert(pars.end(),funcPars.begin(),funcPars.end());
	};
	return pars;
};
//#######################################################################################################################################################
///Gets the constants numbers for wave
std::vector<int> waveset::get_wave_const(
							int 							wave){

	std::vector<int> consts;
	std::vector<int> funcs = get_wave_functions(wave);
	unsigned int nFu = funcs.size();
	for (unsigned int i=0;i<nFu;i++){
		std::vector<int> funcConst = get_function_const(funcs[i]);
		consts.insert(consts.end(),funcConst.begin(),funcConst.end());
	};
	return consts;
};
//########################################################################################################################################################
///Gets all isobars used by wave
std::vector<int> waveset::get_wave_isobars(
							int 							wave){

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
//########################################################################################################################################################
///Gets the isobar paramters for wave
std::vector<int> waveset::get_wave_iso_pars(
							int 							wave){

	std::vector<int> isos = get_wave_isobars(wave);
	std::vector<int> pars;
	for (int i=0;i<isos.size();i++){
		std::vector<int> isopars = get_isobar_pars(isos[i]);
		pars.insert(pars.end(),isopars.begin(),isopars.end());
	};
	return pars;
};
//########################################################################################################################################################
///Gets the isobar constants for wave
std::vector<int> waveset::get_wave_iso_const(
							int 							wave){

	std::vector<int> isos = get_wave_isobars(wave);
	std::vector<int> consts;
	for (int i=0;i<isos.size();i++){
		std::vector<int> isoconst = get_isobar_const(isos[i]);
		consts.insert(consts.end(),isoconst.begin(),isoconst.end());
	};
	return consts;
};
//########################################################################################################################################################
/// Simple getter for function names
std::string waveset::getFunctionName(
							int 							i){

	return _funcNames[i];
};
//########################################################################################################################################################
///Gets the numbers of parameters for each func
std::vector<int> waveset::get_nPars(){

	std::vector<int> nPars;
	nPars.push_back(_borders_par[0]);
	unsigned int i_max = _borders_par.size();
	for (unsigned int i=1;i<i_max;i++){
		int diff = _borders_par[i]-_borders_par[i-1];
		nPars.push_back(diff);
	};
	return nPars;
};
//########################################################################################################################################################
///Gets the numbers of constants for each func
std::vector<int> waveset::get_nConst(){

	std::vector<int> nConst;
	nConst.push_back(_borders_const[0]);
	int i_max = _borders_const.size();
	for (unsigned int i=1;i<i_max;i++){
		int diff = _borders_const[i]-_borders_const[i-1];
		nConst.push_back(diff);
	};
	return nConst;
};
//########################################################################################################################################################
///Gets the paramters numbers for func
std::vector<int> waveset::get_function_pars(
							int 							func){

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
//########################################################################################################################################################
///Get the constant numbers for func
std::vector<int> waveset::get_function_const(
							int 							func){

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
//########################################################################################################################################################
///Get the waves that employ the function number 'func'
std::vector<int> waveset::get_function_waves(
							int 							func){

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
//########################################################################################################################################################
///Simple getter for isobar names
std::string waveset::getIsobarName(
							int 							i){

	return _iso_funcNames[i];
};
//########################################################################################################################################################
///Gets the numbers of paramters for each isobar
std::vector<int> waveset::get_nParsIso(){

	std::vector<int> nPars;
	nPars.push_back(_iso_borders_par[0]);
	unsigned int i_max = _iso_borders_par.size();
	for (unsigned int i=1;i<i_max;i++){
		int diff = _iso_borders_par[i]-_iso_borders_par[i-1];
		nPars.push_back(diff);
	};
	return nPars;
};
//########################################################################################################################################################
///Gets the numbers of constants for each isobar
std::vector<int> waveset::get_nConstIso(){

	std::vector<int> nConst;
	nConst.push_back(_iso_borders_const[0]);
	int i_max = _iso_borders_const.size();
	for (unsigned int i=1;i<i_max;i++){
		int diff = _iso_borders_const[i]-_iso_borders_const[i-1];
		nConst.push_back(diff);
	};
	return nConst;
};
//########################################################################################################################################################
///Gets the numbers of paramters for the 'func' isobar paramterization
std::vector<int> waveset::get_isobar_pars(
							int 							func){

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
//########################################################################################################################################################
///Gets the numbers of constants for the 'func' isobar paramterization
std::vector<int> waveset::get_isobar_const(
							int 							func){

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
//########################################################################################################################################################
///Get the waves that employ a certain isobar
std::vector<int> waveset::get_isobar_waves(
							int 							iso){

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
//########################################################################################################################################################
///Simple getter for parameter names
std::string waveset::getParameterName(
							int 							i){ 	// # of paramter

	return _parNames[i];
};
//########################################################################################################################################################
///Simple getter for constant names
std::string waveset::getConstantName(
							int 							i){ 	// # of constant

	return _constNames[i];
};
//########################################################################################################################################################
///Simple getter for isobar parameter names
std::string waveset::getIsoParName(
							int 							i){ 	// # of isobar paramter

	return  _iso_parNames[i];
};
//########################################################################################################################################################
///Simple getter for isobar constant names
std::string waveset::getIsoConstName(
							int 							i){ 	// # of isobar constant
	return _iso_constNames[i];
};
//########################################################################################################################################################
///Gets the mass bin for a certain m3Pi
int waveset::get_bin(
							double 							mass){

	for(int i=0;i<_nBins;i++){
		if (_binning[i] <= mass and _binning[i+1]>mass){
			return i;
		};
	};
	std::cerr<<"Error: waveset.cxx: get_bin(): Mass not in range"<<std::endl;
	return 0;
};
//########################################################################################################################################################
///Get the first function for each branching-set
std::vector<int> waveset::getFirstBranch(){

	std::vector<int> firstBranch;
	int maxBranch = -1;
	for (int i=0;i<_coupled.size();i++){
		if (_coupled[i] > maxBranch){
			maxBranch = _n_branch[i];
		};
	};
	for (int i=0;i<maxBranch+1;i++){
		for (int j=0;j<_coupled.size();j++){
			if (_coupled[j] == i){
				firstBranch.push_back(_n_branch[j]);
				break;
			};
		};
	};
	return firstBranch;
};
//########################################################################################################################################################
///Updates the number of function-to-wave couplings
void waveset::updateNftw(){

	_nFtw = _funcs_to_waves.size();
};
//########################################################################################################################################################
///Updates the number of employed points
void waveset::updateNpoints(){

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
	};
	_nPoints = nnn;
	_point_to_wave = point_to_wave;
};
//########################################################################################################################################################
///Updates the mass limits for the funtions from the mass limits for the waves
void waveset::updateFuncLims(){

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
//########################################################################################################################################################
//Updates the function spin according to their corresponding wave spin (if there is a mistake (different spins for the same function), nothing will happen)
void waveset::updateFuncSpin(){

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
//########################################################################################################################################################
///Updates the internal definitions for the isobar paramterizations
void waveset::updateIsobar(){

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
//########################################################################################################################################################
///Updates the internal ranges for evaluation
void waveset::update_min_max_bin(){

	_minBin = 0;
//	std::cout<<"Call update_min_max_bin()"<<std::endl;
	_maxBin = _binning.size();
	if (_lowerLims.size() == 0 or _upperLims.size() == 0){
		std::cout<< "Warning: Wave limits not set, omitting internal setting of _minBin and _maxBin"<<std::endl;
		std::cout<< "Warning: Setting _minBin = 0; _maxBin = _binning.size() = "<<_binning.size()<<std::endl;
		_minBin =0;
		_maxBin =_binning.size();
		return;
	};
	double ancMin = _lowerLims[0];
	double ancMax = _upperLims[0];
	for (unsigned int i=0;i<_lowerLims.size();i++){
		if(_lowerLims[i] < ancMin){
			ancMin = _lowerLims[i];
		};
		if(_upperLims[i] > ancMax){
			ancMax = _upperLims[i];
		};
	};
//	std::cout<<ancMin<<"-"<<ancMax<<std::endl;
	if (_binning.size() == 0){
		std::cout<< "Warning: No binning set, cannot determine _minBin, _maxBin" <<std::endl;
	};
	for (unsigned int i=0;i<_binning.size()-1;i++){
		double up = _binning[i+1];
		double low= _binning[i];
		if (ancMin >= low and ancMin <up){
			_minBin = i;
		};
		if (ancMax > low and ancMax <=up){
			_maxBin = i+1;
		};
	};
//	std::cout<<"_minBin: "<<_minBin<<std::endl;
//	std::cout<<"_maxBin: "<<_maxBin<<std::endl;
};
//########################################################################################################################################################
///Handles the branchings after a function/wave coupling was added
void waveset::handle_branchings(
							int 							wave,
							int 							func){

	int border = _borders_waves[wave]-1;
	std::vector<int> coupled;
	std::vector<int> n_branch;
	unsigned int nRel = _coupled.size();
	for (unsigned int i=0; i<nRel;i++){
		coupled.push_back(_coupled[i]);
		n_branch.push_back(_n_branch[i]);
		if (i==border-1){
			n_branch.push_back(-1);
			coupled.push_back(-1);
		};
	};
	if (coupled.size() ==0){
		n_branch.push_back(-1);
		coupled.push_back(-1);
	};
	_coupled=coupled;
	_n_branch = n_branch;
	update_n_cpls();
	update_n_branch();
};
//########################################################################################################################################################
///Updates the number of couplings
void waveset::update_n_cpls(){

	std::vector<int> new_cpl;
	int max_cpl = 0;
	std::vector<int> brs;
	std::vector<int> their_cpls;
	for (int i=i;i<_coupled.size();i++){
		if (_coupled[i] ==-1){
			new_cpl.push_back(max_cpl);
			max_cpl++;
		}else{
			bool already = false;
			for (int j=0;j<brs.size();j++){
				if(brs[j] == _coupled[i]){
					new_cpl.push_back(their_cpls[j]);
					already = true;
					break;
				};
			};
			if (not already){
				new_cpl.push_back(max_cpl);
				brs.push_back(_coupled[i]);
				their_cpls.push_back(max_cpl);
				max_cpl++;
			};
		};
	};
	_nBrCpl = max_cpl;
	_n_cpls = new_cpl;
};
//########################################################################################################################################################
///Updates the branchings
void waveset::update_n_branch(){

	std::vector<int> n_branch;
	int act_branch = 0;
	for (int i=0;i<_coupled.size();i++){
		if (_coupled[i] == -1){
			n_branch.push_back(-1);
		}else{
			n_branch.push_back(act_branch);
			act_branch++;
		};
	};
	_n_branch = n_branch;
	_nBranch = act_branch;
};
//########################################################################################################################################################
///Sets the value for t' fot the corresponding t' bin [tbin] for each constant that is t'
void waveset::updateTprime(
							int 							tbin){

	double tt = _t_binning[tbin] * 0.7 + _t_binning[tbin+1] * 0.3;
//	std::cout<<"set t to: "<<tt<<std::endl;
	for (unsigned int i = 0;i<_const_is_t.size();i++){
		_const[_const_is_t[i]] = tt;
	};
};
//########################################################################################################################################################
///Gives the calss name 'waveset'
std::string className(){

	return "waveset";
};
//########################################################################################################################################################
///Does some internal consistency checks.
bool waveset::checkConsistency(){

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
	if (nErr>0){
		return false;
	};
	return true;
};
//########################################################################################################################################################
///Prints the internal status of the class instance
void waveset::printStatus(){

	std::cout <<"Internal status of waveset:"<<std::endl<<std::endl;
	std::cout<<"WAVES:"<<std::endl;
	std::cout << "_nWaves: " <<_nWaves<<std::endl;
	std::cout << std::endl;
	std::cout << "_waveNames" << std::endl;
	print_vector(_waveNames);
	std::cout<<std::endl<<"_nPoints: "<<_nPoints<<std::endl;
	std::cout << std::endl;
	std::cout<<"_nFtw (number of funtions to waves): "<<_nFtw<<std::endl;
	std::cout<<std::endl;
	std::cout << "_upperLims" << std::endl;
	print_vector(_upperLims);
	std::cout << "_lowerLims" << std::endl;
	print_vector(_lowerLims);
	std::cout << std::endl;
	std::cout << "_L" << std::endl;
	print_vector(_L);
	std::cout<<std::endl<<"_L_iso"<<std::endl;
	print_vector(_L_iso);
	std::cout<<std::endl;
	std::cout<<std::endl<<"_wave_binning_pts"<<std::endl;
	print_vector(_wave_binning_pts);
	std::cout<<std::endl;
	std::cout<<std::endl<<"_point_to_wave"<<std::endl;
	print_vector(_point_to_wave);
	std::cout<<std::endl;
	std::cout<<std::endl<<"_wave_n_binning"<<std::endl;
	print_vector(_wave_n_binning);
	std::cout<<std::endl<<std::endl<<"FUNCTIONS: "<<std::endl;
	std::cout << "_funcs_to_waves" << std::endl;
	print_vector(_funcs_to_waves);
	std::cout<<std::endl<<"_n_cpls"<<std::endl;
	print_vector(_n_cpls);
	std::cout << std::endl;
	std::cout<<"_nPar: "<<_nPar<<std::endl;
	std::cout << std::endl;
	std::cout << "_nFuncs: " << _nFuncs<<std::endl;
	std::cout << std::endl;
	std::cout << "_funcNames" << std::endl;
	print_vector(_funcNames);
	std::cout<<std::endl;
	std::cout << "_funcs"<<std::endl;
	print_vector(_funcs);
	std::cout << std::endl;
	std::cout << "_borders_waves" << std::endl;
	print_vector(_borders_waves);
	std::cout<<std::endl;
	std::cout<<"_maxNpars: "<<_maxNpars<<std::endl;
	std::cout<<std::endl;
	std::cout << "_funcUpperLims" << std::endl;
	print_vector(_funcUpperLims);
	std::cout << "_funcLowerLims" << std::endl;
	print_vector(_funcLowerLims);
	std::cout << std::endl;
	std::cout<<"_L_func"<<std::endl;
	print_vector(_L_func);
	std::cout<<std::endl<<std::endl<<"ISOBARS: "<<std::endl;
	std::cout<<std::endl<<"_iso_to_waves"<<std::endl;
	print_vector(_iso_to_waves);
	std::cout << std::endl;
	std::cout<<std::endl<<"_nIso: "<<_nIso<<std::endl;
	std::cout << std::endl;
	std::cout<<std::endl<<"_iso_funcNames"<<std::endl;
	print_vector(_iso_funcNames);
	std::cout << std::endl;
	std::cout<<std::endl<<"_isos"<<std::endl;
	print_vector(_isos);
	std::cout << std::endl;
	std::cout<<std::endl<<"_maxNparsIso: "<<_maxNparsIso<<"; _maxBinIso: "<<_maxBinIso<<std::endl;
	std::cout << std::endl;
	std::cout<<std::endl<<"_iso_n_binning"<<std::endl;
	print_vector(_iso_n_binning);
	std::cout<<std::endl<<"_iso_binnings"<<std::endl;
	for (unsigned int i=0;i<_iso_binnings.size();i++){
		print_vector(_iso_binnings[i]);
	};
	std::cout << std::endl;
	std::cout<<std::endl<<"_point_borders_wave"<<std::endl;
	print_vector(_point_borders_wave);
	std::cout << std::endl;
	std::cout<<std::endl<<"_L_iso_func"<<std::endl;
	print_vector(_L_iso_func);
	std::cout << std::endl;
	std::cout<<std::endl<<"_L_iso_func"<<std::endl;
	print_vector(_L_iso_func);
	std::cout<<std::endl<<std::endl<<"PARAMETERS & CONSTANTS: "<<std::endl;
	std::cout << "_parNames" << std::endl;
	print_vector(_parNames);
	std::cout << std::endl;
	std::cout << "_borders_par" << std::endl;
	print_vector(_borders_par);
	std::cout << std::endl;
	std::cout << "_constNames" << std::endl;
	print_vector(_constNames);
	std::cout << std::endl;
	std::cout << "_borders_const" << std::endl;
	print_vector(_borders_const);
	std::cout << std::endl;
	std::cout << "_const" << std::endl;
	print_vector(_const);
	std::cout << std::endl;
	std::cout<<std::endl<<"_const_is_t"<<std::endl;
	print_vector(_const_is_t);
	std::cout<<std::endl<<"_iso_parNames"<<std::endl;
	print_vector(_iso_parNames);
	std::cout << std::endl;
	std::cout<<std::endl<<"_iso_borders_par"<<std::endl;
	print_vector(_iso_borders_par);
	std::cout << std::endl;
	std::cout<<std::endl<<"_iso_constNames"<<std::endl;
	print_vector(_iso_constNames);
	std::cout<<std::endl;
	std::cout<<std::endl<<"_iso_borders_const"<<std::endl;
	print_vector(_iso_borders_const);
	std::cout<<std::endl;
	std::cout<<std::endl<<"_iso_const"<<std::endl;
	print_vector(_iso_const);
	std::cout<<std::endl<<std::endl<<"PHASE SPACE: "<<std::endl;
	std::cout<<"_globalPs: "<<_globalPs<<std::endl;
	std::cout << std::endl;
	std::cout << "_wavePs" << std::endl;
	print_vector(_wavePs);
	std::cout << std::endl;
	std::cout<<std::endl<<std::endl<<"BINNING: "<<std::endl;
	std::cout<<"_nBins: "<<_nBins<<std::endl<<std::endl<<"_binning"<<std::endl;
	print_vector(_binning);
	std::cout<<std::endl<<"_minBin: "<<_minBin<<"; _maxBin: "<<_maxBin<<std::endl;
	std::cout<<std::endl<<"_mMin: "<<_mMin<<"; _mMax: "<<_mMax<<std::endl;
	std::cout<<"_nTbin: "<<_nTbin<<std::endl<<std::endl<<"_t_binning"<<std::endl;
	print_vector(_t_binning);
	std::cout<<std::endl;
	print_vector(_eval_tbin);
	std::cout<<std::endl<<std::endl<<"BRANCHING: "<<std::endl;
	std::cout<<"_nBranch: "<<_nBranch<<"; _nBrCpl: "<<_nBrCpl<<std::endl<<"_coupled"<<std::endl;
	print_vector(_coupled);
	std::cout<<std::endl<<"_n_branch"<<std::endl;
	print_vector(_n_branch);
	std::cout<<std::endl<<std::endl<<"INTERNAL: "<<std::endl;
	std::cout<<"_write_out: "<<_write_out<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl<<"_has_isobars: "<<_has_isobars<<std::endl;
};
//########################################################################################################################################################
///Prints the paramters in a 'nice' way
void waveset::printParameters(){

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
//########################################################################################################################################################
///Initializes the text output. Opens a file with name
void waveset::open_output(
							std::string 							name){

	_outStream = new std::ofstream;
	_outStream->open(name.c_str());
	_write_out = true;
};
//########################################################################################################################################################
///Closes the text output
void waveset::close_output(){

	_outStream->close();
	_write_out = false;
};
//########################################################################################################################################################