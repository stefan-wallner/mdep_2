#include"anchor_t.h"
#include<vector>
#include<complex>
#include<fstream>
#include<iostream>
#include<string>

#include<Eigen/Dense>

#include <adolc/adolc.h>  

#include"invert33.h"

anchor_t::anchor_t(): chi2_2d(), _verbose(false),_is_ampl(false){};

template<typename xdouble>
xdouble anchor_t::EvalCP(std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par){ // Evaluates chi2 w/o branchings
	xdouble chi2 = 0.;
	std::vector<std::complex<xdouble> > actCpl = std::vector<std::complex<xdouble> >(_nFtw);
	for (int tbin=0;tbin<_nTbin;tbin++){
		if(_eval_tbin[tbin]){
			updateTprime(tbin);
			for (int i=0;i<_nFtw;i++){
				actCpl[i] = cpl[i+tbin*_nFtw];
			};
			chi2+=EvalTbin(tbin,actCpl,par);
	//		std::cout<<tbin<<":E:"<<EvalTbin(tbin,actCpl,par)<<std::endl;
		};
	};
	return chi2;
};
template double anchor_t::EvalCP(std::vector<std::complex<double> > &cpl,std::vector<double> &par);


template<typename xdouble>
xdouble anchor_t::EvalBranch(std::vector<std::complex<xdouble> >&branch, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par){ // Evaluates chi2 with branchings
	std::vector<std::complex<xdouble> > cplin = std::vector<std::complex<xdouble> >(_nFtw*_nTbin);
	for (int tbin=0;tbin<_nTbin;tbin++){
		if(_eval_tbin[tbin]){
			updateTprime(tbin);
			for (int i =0;i<_nFtw;i++){
				if (_n_branch[i]==-1){
					cplin[i+tbin*_nFtw] = cpl[_n_cpls[i]+tbin*_nBrCpl];
				}else{
					cplin[i+tbin*_nFtw] = cpl[_n_cpls[i]+tbin*_nBrCpl] * branch[_n_branch[i]];
				};
			};
		};
	};
	return EvalCP(cplin,par);
};
template double anchor_t::EvalBranch(std::vector<std::complex<double> >&branch, std::vector<std::complex<double> > &cpl, std::vector<double> &par);

double anchor_t::EvalAutoCpl(std::vector<std::complex<double> > &cpl,std::vector<double> &par){ // gets Chi2 with automatic non anchor couplings (no branchings)
	double chi2 = 0.;
	int nNon = _borders_waves[0];
	std::vector<std::complex<double> > actCpl = std::vector<std::complex<double> >(nNon);
	for (int tbin=0;tbin<_nTbin;tbin++){
		if(_eval_tbin[tbin]){
			updateTprime(tbin);
			for(int i=0;i<nNon;i++){
				actCpl[i] = cpl[i+tbin*nNon];
			};
			chi2+=EvalAutoCplTbin(tbin,actCpl,par);
	//		std::cout << tbin<<":A:"<< EvalAutoCplTbin(tbin,actCpl,par)<<std::endl;
		};
	};
	return chi2;
};

double anchor_t::EvalAutoCplBranch(std::vector<std::complex<double> >&bra, std::vector<std::complex<double> >&cpl, std::vector<double> &par){  // Gets Chi2 for automatically calculated couplings with branchings (the BEST)
	double chi2=0.;
//	std::cout<<"par[1]: "<<par[1]<<std::endl;
	int nNon = _nBrCplAnc;
	std::vector<std::complex<double> > actCpl = std::vector<std::complex<double> >(nNon);
	for(int tbin=0;tbin<_nTbin;tbin++){
		if (_eval_tbin[tbin]){
			updateTprime(tbin);
			for(int i=0;i<nNon;i++){
				actCpl[i] = cpl[i+tbin*nNon];
			};
			std::vector<std::complex<double> > bestcpl = getMinimumCplBra(tbin,bra,actCpl,par);
		//	print_vector(bestcpl);
			std::vector<std::complex<double> > best_cpl_std = std::vector<std::complex<double> >(_nFtw);
			for (int i=0;i<_nFtw;i++){
				if(-1==_n_branch[i]){
					best_cpl_std[i] = bestcpl[_n_cpls[i]];
				}else{
					best_cpl_std[i] = bestcpl[_n_cpls[i]] * bra[_n_branch[i]]; 
				};
			};
			chi2+=EvalTbin(tbin,best_cpl_std,par);
		};
	};
	return chi2;
};

template<typename xdouble>
xdouble anchor_t::EvalTbin(int tbin, std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par){ // gets the chi2 for a single t' bin
	updateTprime(tbin);
	xdouble chi2 = 0.;
	for (int bin=_minBin; bin<_maxBin; bin++){
		chi2+=EvalBin(tbin,bin,cpl,par);
//		std::cout << bin<<"   "<< EvalBin(tbin,bin,cpl,par) <<std::endl;
		if(_verbose){      
			std::cout << "────────────────────────────────────────────────" <<std::endl;
		};
	};
	return chi2;
};
template double anchor_t::EvalTbin(int tbin, std::vector<std::complex<double> > &cpl,std::vector<double> &par);

double anchor_t::EvalAutoCplTbin(int tbin, std::vector<std::complex<double> > &cpl, std::vector<double> &par){ // Gets the chi2 for a certain t' bin, with automatically calculated couplings
	std::vector<std::complex<double> > mincpl = getMinimumCpl(tbin,cpl,par);
	return EvalTbin(tbin,mincpl,par);
};

template<typename xdouble>
xdouble anchor_t::EvalBin(int tbin,int bin,std::vector<std::complex<xdouble> >&cpl,std::vector<xdouble> &par){ // Gets the Chi2 for a single t' and m3pi bin
	double mass = (_binning[bin] + _binning[bin+1])/2; // Eval at bin center.
	std::vector<xdouble> deltas = delta(tbin,bin,mass, cpl,par);
	xdouble chi2 = 0.;
	for (int i=0;i<2*_nWaves-1;i++){
		int iWave = (i+1)/2;
		if (mass >= _lowerLims[iWave] and mass < _upperLims[iWave]){
			if(_verbose){
				std::cout <<"i: "<<i<<" j: "<<i<<" bin: "<<bin<< " mass: "<<mass<<" chi2+: "<<2.*deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<std::endl;
			};
//			std::cout<<"le_addite:D"<<pow(deltas[i],2.)*_coma[tbin][bin][i][i]<<std::endl;
			chi2+= deltas[i]*deltas[i]*_coma[tbin][bin][i][i];
			for (int j=0;j<i;j++){
				int jWave = (j+1)/2;
				if(mass >= _lowerLims[jWave] and mass < _upperLims[jWave]){
					if (_verbose){
						std::cout <<"i: "<<i<<" j: "<<j<<" bin: "<<bin<< " mass: "<<mass<<" chi2+: "<<2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]<<std::endl;
					};
//					std::cout<<"le_addite: "<<2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]<<std::endl;
					chi2+=2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]; // Factor 2. because _coma[][][][] is symmetric
				};	
			};
		};
	};
	return chi2;
};
template double anchor_t::EvalBin(int tbin,int bin,std::vector<std::complex<double> >&cpl,std::vector<double> &par);

template<typename xdouble>
std::vector<xdouble> anchor_t::delta(int tbin, int bin,double mass, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par){ // Returns f(m,...) - data[...] for each SDM entry in the fit
	std::vector<std::vector<std::complex<xdouble> > > iso_eval;
	std::vector<std::complex<xdouble> > ampls = amps(mass, cpl, par, iso_eval);
	std::vector<xdouble> del = std::vector<xdouble> (2*_nWaves-1);
	std::complex<xdouble> divider = std::complex<xdouble>(1.,0.); // When amplitudes are fitted, this is set to |ampls[0]|
	if(_is_ampl){
		divider = std::complex<xdouble>(pow(std::norm(ampls[0]),.5),0.);
	};
	del[0]=std::norm(ampls[0]/divider) - _data[tbin][bin][0];
	if(_verbose){
		std::cout<<"anc: "<<std::norm(ampls[0])<<" - "<<_data[tbin][bin][0]<<std::endl;
	};
	for (int i = 1; i<_nWaves;i++){
		std::complex<xdouble> inter = ampls[0]*std::conj(ampls[i])/divider;
		if(_verbose){
			std::cout<<"real"<<i<<": "<<real(inter)<<" - "<<_data[tbin][bin][2*i-1]<<std::endl;
			std::cout<<"imag"<<i<<": "<<imag(inter)<<" - "<<_data[tbin][bin][2*i]<<std::endl;
		};
		del[2*i-1]=real(inter) - _data[tbin][bin][2*i-1]; // real part
		del[2*i]=imag(inter) - _data[tbin][bin][2*i];    // imag part
	};
	return del;
};
template std::vector<double> anchor_t::delta(int tbin, int bin,double mass, std::vector<std::complex<double> > &cpl, std::vector<double> &par);

std::vector<std::vector<double> > anchor_t::getPlots(int tbin, std::vector<std::complex<double> >&cpl, std::vector<double> &par){ ///// DOES NOT WORK !!!!!
	updateTprime(tbin);
	std::vector<std::vector<double> > plots = std::vector<std::vector<double> >(2*_nWaves -1); // Store everything for one data set in one vector, with the structure [fit, data, error, fit, data, error, fit, ... ]
	for(int bin=0;bin<_nBins;bin++){
		std::vector<std::vector<std::complex<double> > > iso_eval; // If isobars were used
		double mass = (_binning[bin]+_binning[bin+1])/2.;
		std::vector<std::complex<double> > ampls = amps(mass, cpl ,par,iso_eval);
		std::complex<double> divider = std::complex<double>(1.,0.); // When amplitudes are fitted, this is set to |ampls[0]|
		if(_is_ampl){
			divider = std::complex<double>(pow(std::norm(ampls[0]),.5),0.);
		};
		plots[0].push_back(std::norm(ampls[0]/divider));
		plots[0].push_back(_data[tbin][bin][0]);
		plots[0].push_back(_coma[tbin][bin][0][0]); // Use simple diagonal errors here.
		for (int wave=1;wave<_nWaves;wave++){
			std::complex<double> inter = ampls[0]*std::conj(ampls[wave])/divider;
			plots[2*wave-1].push_back(real(inter));
			plots[2*wave-1].push_back(_data[tbin][bin][2*wave-1]);
			plots[2*wave-1].push_back(_coma[tbin][bin][2*wave-1][2*wave-1]);
			plots[2*wave  ].push_back(real(inter));
			plots[2*wave  ].push_back(_data[tbin][bin][2*wave  ]);
			plots[2*wave  ].push_back(_coma[tbin][bin][2*wave  ][2*wave  ]);

		};
	};
	return plots;
};

#ifdef ADOL_ON // Enable adouble operations, if the auto diff package is loaded
template adouble anchor_t::EvalCP(std::vector<std::complex<adouble> > &cpl,std::vector<adouble> &par);
template adouble anchor_t::EvalBranch(std::vector<std::complex<adouble> >&branch, std::vector<std::complex<adouble> > &cpl, std::vector<adouble> &par);
template adouble anchor_t::EvalTbin(int tbin, std::vector<std::complex<adouble> > &cpl,std::vector<adouble> &par);
template adouble anchor_t::EvalBin(int tbin,int bin,std::vector<std::complex<adouble> >&cpl,std::vector<adouble> &par);
template std::vector<adouble> anchor_t::delta(int tbin, int bin,double mass, std::vector<std::complex<adouble> > &cpl, std::vector<adouble> &par);
#endif//ADOL_ON

void anchor_t::setBinning(std::vector<double> binning){ // sets the binning in m3pi
	_binning = binning;
	_nBins = binning.size()-1;
	_mMin = binning[0];
	_mMax = binning[_nBins];
	int minBin = 0;
	int maxBin = _nBins-1;
	double ancMin = _lowerLims[0];
	double ancMax = _upperLims[0];
	bool upToSet = true;
	bool lowToSet= true;
	for (int bin =0; bin<_nBins;bin++){
		double bug = binning[bin];	//  Binobergrenze
		double bog = binning[bin+1];	//  BinUntergrenze
		double bc = (bog+bug)/2;
		if (bc >= ancMin and lowToSet){
			minBin = bin;
			lowToSet=false;
		};
		if (bog > ancMax and upToSet){
			maxBin = bin-1;
			upToSet=false;
		};
	};
	_minBin = minBin;
	_maxBin = maxBin;
};

std::string anchor_t::className(){
	return "anchor_t";
};

void anchor_t::setTbinning(std::vector<double> binning){ // sets the binning in t'
	_t_binning = binning;
	_nTbin = binning.size() - 1;
	_data = std::vector<std::vector<std::vector<double> > >(_nTbin);
	_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_nTbin);
	if (_eval_tbin.size() != _nTbin){
		_eval_tbin = std::vector<bool>(_nTbin,true);
	};
};

int anchor_t::get_bin(double mass){ // gets the mass bin for a certain m3Pi
	for(int i=0;i<_nBins;i++){
		if (_binning[i] <= mass and _binning[i+1]>mass){
			return i;
		};
	};
	std::cerr<<"Error: anchor_t.cxx: get_bin(): Mass not in range"<<std::endl;
	return 0;
};
std::vector<double> anchor_t::get_data(int tbin,int bin){ // returns the data for a t' bin
	return _data[tbin][bin];
};
std::vector<std::vector<double> > anchor_t::get_coma(int tbin,int bin){ // returns the coma for a t' bin
	return _coma[tbin][bin];
};

void anchor_t::loadData(int tbin, const char* dataFile){ // Loads the data for a t' bin
	_data[tbin] = std::vector<std::vector<double> >();
	std::fstream data(dataFile,std::ios_base::in);
	double val;
	unsigned int nDat = 2*_nWaves-1;
	std::vector<double> data_bin;
	while(data >> val){
		data_bin.push_back(val);
		if (data_bin.size() == nDat){
			_data[tbin].push_back(data_bin);
			data_bin = std::vector<double>();
		};		
	};
	if (_nBins != _data[tbin].size()){
		std::cout << "Warning: _nBins != _data.size()"<<std::endl;
	}else{		
		std::cout << "File delivered the right size for _data"<<std::endl;	
	};
};

void anchor_t::loadComa(int tbin, const char* comaFile){ // Loads the (inverted) covariance matrix for a t' bin from 'comaFile'
	_coma[tbin]=std::vector<std::vector<std::vector<double> > >();
	std::fstream data(comaFile,std::ios_base::in);
	double val;
	unsigned int nDat = 2*_nWaves-1;
	std::vector<std::vector<double> > coma_bin;
	std::vector<double> coma_line;
	while(data>>val){
		coma_line.push_back(val);
		if (coma_line.size() == nDat){
			coma_bin.push_back(coma_line);
			coma_line = std::vector<double>();
		};
		if (coma_bin.size() == nDat){
			_coma[tbin].push_back(coma_bin);
			coma_bin = std::vector<std::vector<double> >();
		};
	};
	if (_nBins != _coma[tbin].size()){
		std::cout << "Warning: _nBins != _coma.size()" << std::endl;
	}else{
		std::cout<< "File delivered the right size for _coma"<<std::endl;
	};	
};

void anchor_t::setEvalTbin(int i, bool flag){ // Switches single t' bins on/off 
	_eval_tbin[i] = flag;
};

void anchor_t::nullify(){ // Sets all data to 0.
	for (int tbin = 0; tbin < _nTbin;tbin++){
		for (int bin =0; bin<_nBins;bin++){
			int nDat = 2*_nWaves-1;
			for (int i=0;i<nDat;i++){
				_data[tbin][bin][i]=0.;	
			};	
		};
	};
};

void anchor_t::add_func(int i, bool is_t_dep){ // Adds a function
	chi2_2d::add_func(i);
	if (is_t_dep){
		_const_is_t.push_back(_borders_const[_borders_const.size()-1]-1);
	};
};

void anchor_t::setConst(int i,double con){ // Sets a constant
	chi2_2d::setConst(i,con);
	for (unsigned int j=0;j<_const_is_t.size();j++){
		if (i == _const_is_t[j]){
			std::cout<<"Warning: Trying to set _const["<<i<<"] which is t' and will be overwritten"<<std::endl;
		};
	};
};
void anchor_t::add_func_to_wave(int wave, int func){ // Adds a function paramtetrization to a wave amplitude
	chi2_2d::add_func_to_wave(wave,func);
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
void anchor_t::couple_funcs(int i1,int i2){ // Couples two functions: {C1(t), C2(t)} -> {B1*C(t), B2*C(t)} 
	if (i1 == i2){
		std::cerr<<"Error: Trying to couple two times the same functions"<<std::endl;
		return;
	};
	int branch1 = _coupled[i1];
	int branch2 = _coupled[i2];
	if (branch1 == branch2 and branch1 != -1){
		return;
	};
	if (branch1 == -1 and branch2 == -1){
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
	}else if(branch1 == -1){
		_coupled[i1] = branch2;
	}else if(branch2 == -1){
		_coupled[i2] = branch1;
	}else{
		std::cerr<<"function["<<i1<<"] and function["<<i2<<"] are already in different branchings"<<std::endl;
	};
	update_n_cpls();
	update_n_branch();
};

void anchor_t::setVerbose(bool in){
	_verbose = in;
};
/*
		┌───────────────────────────────────────┐
		│	Es folgt eine Indexschlacht	│
		└───────────────────────────────────────┘
*/
AandB anchor_t::get_AB(int tbin,std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par){// Gets the coefficients A and B for Chi2 = x^T*A*x + B*x + C
	int nCplAnc = _borders_waves[0]; // Number of couplings for the anchor wave
	int nNon = _nFtw - nCplAnc;
	AandB AB(2*nNon);
	std::vector<double> lines = std::vector<double>(2*_nWaves-2,0.);
	for(int bin=_minBin;bin<_maxBin;bin++){
		double m = (_binning[bin]+_binning[bin+1])/2;
		std::vector<std::complex<double> > func = funcs(m,par);


		std::vector<double> phase = phase_space(m);
		std::complex<double> ampAnc = std::complex<double>(0.,0.);
		for(int i=0;i<nCplAnc;i++){
			ampAnc+= anchor_cpl[i]*func[_funcs_to_waves[i]]*phase[0];
		};
		if(_is_ampl){ // Divide ampAnc by |ampAnc| to get the amplitude with phase 0. // Need to be tested
			ampAnc/=pow(std::norm(ampAnc),.5);
		};
		double RR = ampAnc.real();
		double II = ampAnc.imag();
		for (int i =0;i<2*_nWaves-2;i++){
			double val = (RR*RR+II*II - _data[tbin][bin][0])*_coma[tbin][bin][i+1][0];
			for (int j=1;j<2*_nWaves-1;j++){
				val+=-_data[tbin][bin][j]*_coma[tbin][bin][i+1][j];
			};
			lines[i]=val;			
		};
		double intAnc = RR*RR+II*II;
		int wave1 = 1;	// Start with 1, leave anchor wave out
		int up1 = _borders_waves[1];
		for (unsigned int nf1 = nCplAnc;nf1<_nFtw;nf1++){
			if (nf1==up1){
				wave1+=1;
				up1 = _borders_waves[wave1];
			};
			int f1 = _funcs_to_waves[nf1];
			int ind1=nf1-nCplAnc;
			int wave2=1;
			int up2=_borders_waves[1];

			std::complex<double> AB1 = ampAnc * conj(func[f1]) * phase[wave1];

			double r1 = func[f1].real()*phase[wave1];
			double i1 = func[f1].imag()*phase[wave1];
			//// (RR + j II)*(R - j I)*(r - j i) = (RR R r - RR I i - II I r + II R i)+j(II R r + II I i - RR r I - RR R i)

			double RR1 = RR*r1+II*i1; // = Coefficient of the real part of the coupling 1 in the real part of wave one
			double RI1 =-RR*i1-II*r1; // = Coefficient of the imag part of the coupling 1 in the real part of wave one
			double IR1 = II*r1-RR*i1; // = Coefficient of the real part of the coupling 1 in the imag part of wave one
			double II1 = II*i1-RR*r1; // = Coefficient of the imag part of the coupling 1 in the imag part of wave one

			AB.B[2*ind1  ] += AB1.real() * lines[2*wave1-2] * 2;
			AB.B[2*ind1  ] += AB1.imag() * lines[2*wave1-1] * 2;

			AB.B[2*ind1+1] +=-AB1.real() * lines[2*wave1-1] * 2;
			AB.B[2*ind1+1] += AB1.imag() * lines[2*wave1-2] * 2;

			for (unsigned int nf2= nCplAnc;nf2<_nFtw;nf2++){
				if(nf2==up2){
					wave2+=1;
					up2 = _borders_waves[wave2];
				};
				int f2 = _funcs_to_waves[nf2];
				int ind2=nf2-nCplAnc;


				std::complex<double> AB2 = ampAnc * conj(func[f2]) * phase[wave2];

				double r2 = func[f2].real()*phase[wave2];
				double i2 = func[f2].imag()*phase[wave2];
				//// (RR + j II)*(R - j I)*(r - j i) = (RR R r - RR I i - II I r + II R i)+j(II R r + II I i - RR r I - RR R i)
				double RR2 = RR*r2+II*i2;
				double RI2 =-RR*i2-II*r2;
				double IR2 = II*r2-RR*i2;
				double II2 = II*i2-RR*r2;
//				std::cout<<"_coma[tbin][bin][2*wave1-1][2*wave2-1]: "<<_coma[tbin][bin][2*wave1-1][2*wave2-1]<<std::endl;
//				std::cout<< "AB1.real(): "<<AB1.real()<<std::endl;
//				std::cout<< "AB2.real(): "<<AB2.real()<<std::endl;

				AB.A(2*ind1  ,2*ind2  )+= AB1.real() * _coma[tbin][bin][2*wave1-1][2*wave2-1] * AB2.real();
				AB.A(2*ind1  ,2*ind2  )+= AB1.real() * _coma[tbin][bin][2*wave1-1][2*wave2  ] * AB2.imag();
				AB.A(2*ind1  ,2*ind2  )+= AB1.imag() * _coma[tbin][bin][2*wave1  ][2*wave2-1] * AB2.real();
				AB.A(2*ind1  ,2*ind2  )+= AB1.imag() * _coma[tbin][bin][2*wave1  ][2*wave2  ] * AB2.imag();

				AB.A(2*ind1  ,2*ind2+1)+=-AB1.real() * _coma[tbin][bin][2*wave1-1][2*wave2  ] * AB2.real();
				AB.A(2*ind1  ,2*ind2+1)+= AB1.real() * _coma[tbin][bin][2*wave1-1][2*wave2-1] * AB2.imag();
				AB.A(2*ind1  ,2*ind2+1)+=-AB1.imag() * _coma[tbin][bin][2*wave1  ][2*wave2  ] * AB2.real();
				AB.A(2*ind1  ,2*ind2+1)+= AB1.imag() * _coma[tbin][bin][2*wave1  ][2*wave2-1] * AB2.imag();

				AB.A(2*ind1+1,2*ind2  )+=-AB1.real() * _coma[tbin][bin][2*wave1  ][2*wave2-1] * AB2.real();
				AB.A(2*ind1+1,2*ind2  )+=-AB1.real() * _coma[tbin][bin][2*wave1  ][2*wave2  ] * AB2.imag();
				AB.A(2*ind1+1,2*ind2  )+= AB1.imag() * _coma[tbin][bin][2*wave1-1][2*wave2-1] * AB2.real();
				AB.A(2*ind1+1,2*ind2  )+= AB1.imag() * _coma[tbin][bin][2*wave1-1][2*wave2  ] * AB2.imag();

				AB.A(2*ind1+1,2*ind2+1)+= AB1.real() * _coma[tbin][bin][2*wave1  ][2*wave2  ] * AB2.real();
				AB.A(2*ind1+1,2*ind2+1)+=-AB1.real() * _coma[tbin][bin][2*wave1  ][2*wave2-1] * AB2.imag();
				AB.A(2*ind1+1,2*ind2+1)+=-AB1.imag() * _coma[tbin][bin][2*wave1-1][2*wave2  ] * AB2.real();
				AB.A(2*ind1+1,2*ind2+1)+= AB1.imag() * _coma[tbin][bin][2*wave1-1][2*wave2-1] * AB2.imag();
			};
		};
	};
	return AB;
};

std::vector<std::complex<double> > anchor_t::getMinimumCpl(int tbin,std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par){	// Gets the bes couplings without branchings
	int nCplAnc = _borders_waves[0]; // Number of couplings for the anchor wave
	int nNon = _nFtw - nCplAnc;
	std::vector<std::complex<double> > cpl = std::vector<std::complex<double> >(_nFtw);
	for (int i=0; i<nCplAnc; i++){
		cpl[i] = anchor_cpl[i];
	};
	AandB AB = get_AB(tbin,anchor_cpl,par);
	Eigen::VectorXd D =  AB.A.ldlt().solve(AB.B); // (*)

	for (int i=0;i<nNon;i++){
		// Short explanation:
		// ┌────────────────────────────────────────────────────────────────────┐
		// │C^T * A * C + B * C + const. = Chi2 is a quadratic form		│
		// │To find the minimum: dChi2/dC_i = 0. for all i in {0,...,2*nNon-1}	│	
		// │=> (A*C)_i + (C^T*A)_i + B_i = 0.					│
		// │Since in the Chi2, the assymmetric part of A cancels out anyway:	│
		// │=> 2*A*C = -B (**)							│
		// │D solves A*D = B (*)						│
		// │D = -2*C								│
		// │=> C solves (**)							│
		// └────────────────────────────────────────────────────────────────────┘
		cpl[nCplAnc+i]=std::complex<double>(-D(2*i)/2.,-D(2*i+1)/2.); 
	};
	return cpl;
};

std::vector<std::complex<double> > anchor_t::getMinimumCplBra(int tbin, std::vector<std::complex<double> > &branch, std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par){ // Calculated all non anchor couplings (no need for fitting)
	if (0==_nBranch){ // Do not do the complicated stuff, when no branchings are used
		return getMinimumCpl(tbin,anchor_cpl,par);
	}else{
		int nCplAnc = _borders_waves[0];
		int nNon = _nFtw - nCplAnc;
		std::vector<std::complex<double> > cplAncBr = std::vector<std::complex<double> >(nCplAnc);
		std::vector<std::complex<double> > ccppll(_nBrCpl);
		for (int i=0;i<_nBrCplAnc;i++){
			ccppll[i] = anchor_cpl[i];
		};
		for (int i=0;i<nCplAnc;i++){
			int ncpl_act = _n_cpls[i];
			int nbranch_act = _n_branch[i];
			if (-1 == nbranch_act){
				cplAncBr[i] = anchor_cpl[ncpl_act];
			}else{
				cplAncBr[i] = anchor_cpl[ncpl_act] * branch[nbranch_act];
			};
		};
		AandB AB = get_AB(tbin,cplAncBr,par);
		AandB ABprime(2*(_nBrCpl - _nBrCplAnc));
		for (int i =0;i<nNon;i++){ // Reshuffle the coefficients 
			int i_tot = i+nCplAnc;
			int i_cpl_tot = _n_cpls[i_tot];
			int i_cpl = i_cpl_tot - _nBrCplAnc;
			int i_br  = _n_branch[i_tot];
	//		std::cout<<"────────────────────────────────────────────────"<<std::endl;
	//		std::cout<<i<<" "<<i_tot<<" "<<i_cpl_tot<<" "<<i_cpl<<" "<<i_br<<std::endl;		
			if (-1 == i_br){ // If the funcion is uncoupled, just take the old coefficients
				ABprime.B(2*i_cpl  )+= AB.B(2*i  );
				ABprime.B(2*i_cpl+1)+= AB.B(2*i+1);
			}else if(i_cpl>-1){ // If the function is not coupled to the anchor wave, add it to the corresponding coupling
				ABprime.B(2*i_cpl  )+= AB.B(2*i  )*branch[i_br].real()   +   AB.B(2*i+1)*branch[i_br].imag();
				ABprime.B(2*i_cpl+1)+=-AB.B(2*i  )*branch[i_br].imag()   +   AB.B(2*i+1)*branch[i_br].real();
			};
			for (int j=0; j<nNon;j++){
				int j_tot = j+nCplAnc;
				int j_cpl_tot = _n_cpls[j_tot];
				int j_cpl = j_cpl_tot -_nBrCplAnc;
				int j_br  = _n_branch[j_tot];
				if ((-1 == i_br or i_cpl>-1)and(-1 == j_br or j_cpl>-1)){
					std::complex<double> ifak(1.,0.);
					std::complex<double> jfak(1.,0.);
					if (i_br >-1){
						ifak = branch[i_br];
					};
					if (j_br >-1){
						jfak = branch[j_br];
					};
					ABprime.A(2*i_cpl  ,2*j_cpl  )+=	 AB.A(2*i  ,2*j  )*ifak.real()*jfak.real() +
										 AB.A(2*i  ,2*j+1)*ifak.real()*jfak.imag() +
										 AB.A(2*i+1,2*j  )*ifak.imag()*jfak.real() +
										 AB.A(2*i+1,2*j+1)*ifak.imag()*jfak.imag() ;
					ABprime.A(2*i_cpl  ,2*j_cpl+1)+=	-AB.A(2*i  ,2*j  )*ifak.real()*jfak.imag() +
										 AB.A(2*i  ,2*j+1)*ifak.real()*jfak.real() +
										-AB.A(2*i+1,2*j  )*ifak.imag()*jfak.imag() +
										 AB.A(2*i+1,2*j+1)*ifak.imag()*jfak.real() ;
					ABprime.A(2*i_cpl+1,2*j_cpl  )+=	-AB.A(2*i  ,2*j  )*ifak.imag()*jfak.real() +
										-AB.A(2*i  ,2*j+1)*ifak.imag()*jfak.imag() +
										 AB.A(2*i+1,2*j  )*ifak.real()*jfak.real() +
										 AB.A(2*i+1,2*j+1)*ifak.real()*jfak.imag() ;
					ABprime.A(2*i_cpl+1,2*j_cpl+1)+=	 AB.A(2*i  ,2*j  )*ifak.imag()*jfak.imag() +
										-AB.A(2*i  ,2*j+1)*ifak.imag()*jfak.real() +
										-AB.A(2*i+1,2*j  )*ifak.real()*jfak.imag() +
										 AB.A(2*i+1,2*j+1)*ifak.real()*jfak.real() ;
	//				std::cout<<ifak<<jfak<<std::endl;
				}else if(-1 == i_br or i_cpl>-1){			
					std::complex<double> ifak(1.,0.);
					if (i_br>-1){
						ifak = branch[i_br];
					};
					std::complex<double> cpl_j = branch[j_br]*anchor_cpl[j_cpl_tot];
					ABprime.B(2*i_cpl  )+=   AB.A(2*i  ,2*j  )*ifak.real()*cpl_j.real() +
								 AB.A(2*i  ,2*j+1)*ifak.real()*cpl_j.imag() +
								 AB.A(2*i+1,2*j  )*ifak.imag()*cpl_j.real() +
								 AB.A(2*i+1,2*j+1)*ifak.imag()*cpl_j.imag() ;
					ABprime.B(2*i_cpl+1)+=	-AB.A(2*i  ,2*j  )*ifak.imag()*cpl_j.real() +
								-AB.A(2*i  ,2*j+1)*ifak.imag()*cpl_j.imag() +
								 AB.A(2*i+1,2*j  )*ifak.real()*cpl_j.real() +
								 AB.A(2*i+1,2*j+1)*ifak.real()*cpl_j.imag() ;
	//				std::cout<<"Nirrr"<<std::endl;
				}else if(-1 == j_br or j_cpl>-1){	
					std::complex<double> jfak(1.,0.);
					if (j_br>-1){
						jfak = branch[j_br];
					};
					std::complex<double> cpl_i = branch[i_br]*anchor_cpl[i_cpl_tot];
					ABprime.B(2*j_cpl  )+=	 AB.A(2*i  ,2*j  )*cpl_i.real()*jfak.real() +
								 AB.A(2*i  ,2*j+1)*cpl_i.real()*jfak.imag() +
								 AB.A(2*i+1,2*j  )*cpl_i.imag()*jfak.real() +
								 AB.A(2*i+1,2*j+1)*cpl_i.imag()*jfak.imag() ;
					ABprime.B(2*j_cpl+1)+=	-AB.A(2*i  ,2*j  )*cpl_i.real()*jfak.imag() +
								 AB.A(2*i  ,2*j+1)*cpl_i.real()*jfak.real() +
								-AB.A(2*i+1,2*j  )*cpl_i.imag()*jfak.imag() +
								 AB.A(2*i+1,2*j+1)*cpl_i.imag()*jfak.real() ;
	//				std::cout<<"Njrrr"<<std::endl;
				};// If both functions are coupled to the anchor wave, they just add to the constant C, (Chi2 = x^TAx + B^Tx +C) and do not change the position of the minimum.
			};
		};
	//	std::vector<double> les_as;
	//	for (int i=0;i<2*(_nBrCpl - _nBrCplAnc);i++){
	//		les_as.push_back(ABprime.A(i,i));
	//	};
	//	print_vector(les_as);

		Eigen::VectorXd D =  ABprime.A.ldlt().solve(ABprime.B); 
		for (int i=_nBrCplAnc;i<_nBrCpl;i++){
			int irel = i-_nBrCplAnc;
			ccppll[i] = std::complex<double>(-D(2*irel)/2,-D(2*irel+1)/2);
		};
	//	std::cout<<"made it through"<<std::endl;
		return ccppll;
	};
};

void anchor_t::updateTprime(int tbin){ // Sets the value for t' fot the corresponding t' bin [tbin] for each constant that is t' 
	double tt = _t_binning[tbin] * 0.7 + _t_binning[tbin+1] * 0.3;
//	std::cout<<"set t to: "<<tt<<std::endl;
	for (unsigned int i = 0;i<_const_is_t.size();i++){
		_const[_const_is_t[i]] = tt;
	};
};

void anchor_t::update_n_cpls(){ // updates the number of couplings
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
	_nBrCplAnc = new_cpl[_borders_waves[0]-1]+1;
	_n_cpls = new_cpl;
};
void anchor_t::update_n_branch(){ // Updates the branchings
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

void anchor_t::printStatus(){ // Prints the internal status
	chi2_2d::printStatus();
	std::cout<<std::endl<< "_nBranch: "<<_nBranch<<std::endl;
	std::cout<<std::endl<< "_nBrCpl: "<<_nBrCpl<<std::endl;
	std::cout<<std::endl<< "_nBrCplAnc: "<<_nBrCplAnc<<std::endl;
	std::cout<<std::endl<< "_coupled:"<<std::endl;
	print_vector(_coupled);
	std::cout<<std::endl<< "_n_cpls:"<<std::endl;	
	print_vector(_n_cpls);
	std::cout<<std::endl<< "_n_branch:"<<std::endl;
	print_vector(_n_branch);
	std::cout<<std::endl<< "_funcs_to_waves:"<<std::endl;
	print_vector(_funcs_to_waves);
	std::cout<<std::endl<<"_binning:"<<std::endl;
	print_vector(_binning);
	std::cout<<std::endl<<"_t_binning:"<<std::endl;
	print_vector(_t_binning);
	std::cout<<std::endl<<"_const_is_t:"<<std::endl;
	print_vector(_const_is_t);
};

int anchor_t::getNtot(){ // total number pf parameters (including non achor couplings)
	return 2*getNcpl() + getNpar() + 2*getNbra() + getNiso();
};
int anchor_t::getNtotAnc(){ // total number of parameters
	return 2*getNanc() + getNpar() + 2*getNbra() + getNiso();
};
int anchor_t::getNcpl(){ // number of total couplings 
	return _nBrCpl * _nTbin;
};
int anchor_t::getNanc(){ // number of total couplings in the anchor wave
	return _nBrCplAnc * _nTbin;
};

int anchor_t::getNpar(){ // returns the number of parameters (without branchings or couplings)
	if (0==_borders_par.size()){
		return 0;
	};
	return _borders_par[_borders_par.size() -1];
};
int anchor_t::getNbra(){ // Number of sets coupled by branchings
	return _nBranch;
};
int anchor_t::getNiso(){ // returns the number of parameters for isobars
	int siz = _iso_borders_par.size();
	if (0==siz){
		return 0;
	};
	return _iso_borders_par[siz-1];
};

int anchor_t::getNtBin(){ // Gives the number of t' bins
	return _nTbin;
};

std::vector<int> anchor_t::getFirstBranch(){// Get the first function for each set 
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

void anchor_t::set_is_ampl(bool is_ampl){ // Set the mode to amplitude fitting, instead of anchor wave
	_is_ampl = is_ampl;
};

std::vector<std::complex<double> > anchor_t::get_branchings(std::vector<std::complex<double> > cpl,std::vector<double> par){ 
	// Gets branchings as follows:
	//	- Take Anchor cpls to get couplings for each t-bin
	//	- Average over the tbins for all couplings
	//	- Determine the branchings from these averages
	std::vector<std::complex<double> > brchs = std::vector<std::complex<double> >(_nBranch);
	std::vector<std::complex<double> > all_cpls = std::vector<std::complex<double> >(_nFtw);
	int active_t_bins=0;
	int cpl_used=0;
	for (int tbin=0;tbin<_nTbin;tbin++){ // Average couplings over all t' bins, after rotating to phase of cpl[0] to 0
		if (_eval_tbin[tbin]){
			std::vector<std::complex<double> > cpl_t = std::vector<std::complex<double> >(_nBrCplAnc);
			for (int i=0;i<_nBrCplAnc;i++){
				cpl_t[i] = cpl[cpl_used];
				cpl_used++;
			};
			active_t_bins++;
			updateTprime(tbin);
			std::vector<std::complex<double> > all_cpls_t= getMinimumCpl(tbin,cpl_t,par);
			std::complex<double> phase = all_cpls_t[0]/abs(all_cpls_t[0]);
			for (int i=0;i<_nFtw;i++){
				all_cpls[i]+=all_cpls_t[i]/phase;
			};
		};
	};
	for (int i=0;i<_nFtw;i++){ // Divide by number of t' bins
		all_cpls[i]/=active_t_bins;
	};
	for (unsigned int j=0;j<_n_branch.size();j++){ // Write coupling-average to corresponding branching
		if (_n_branch[j] >= 0){
			brchs[_n_branch[j]] = all_cpls[j];
		};
	};
	return brchs;
};

