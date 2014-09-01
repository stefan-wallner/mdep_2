#include"anchor_t.h"
#include<vector>
#include<complex>
#include<fstream>
#include<iostream>
#include<string>

#include"cholesky.h"

#include "adolc/adolc.h"

#include"invert33.h"

anchor_t::anchor_t(): waveset(),_is_ampl(false){};

#ifdef ADOL_ON 
std::complex<adouble> operator*(std::complex<adouble> a,double d){// Define std::complex<adouble> * double 
	return std::complex<adouble>(a.real()*d,a.imag()*d);
};
#endif//ADOL_ON

template<typename xdouble>
xdouble anchor_t::EvalCP(std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par, std::vector<xdouble> &iso_par){ // Evaluates chi2 w/o branchings
	xdouble chi2 = 0.;
	std::vector<std::complex<xdouble> > actCpl = std::vector<std::complex<xdouble> >(_nFtw);
	for (int tbin=0;tbin<_nTbin;tbin++){
		if(_eval_tbin[tbin]){
			updateTprime(tbin);
			for (int i=0;i<_nFtw;i++){
				actCpl[i] = cpl[i+tbin*_nFtw];
			};
			chi2+=EvalTbin(tbin,actCpl,par,iso_par);
	//		std::cout<<tbin<<":E:"<<EvalTbin(tbin,actCpl,par)<<std::endl;
//			std::cout<<chi2<<std::endl;
		};
	};
	if (_write_out){
		*_outStream <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double anchor_t::EvalCP(std::vector<std::complex<double> > &cpl,std::vector<double> &par, std::vector<double> &iso_par);


template<typename xdouble>
xdouble anchor_t::EvalBranch(std::vector<std::complex<xdouble> >&branch, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par){ // Evaluates chi2 with branchings
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
	xdouble chi2 = EvalCP(cplin,par, iso_par);
	if (_write_out){
		*_outStream <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double anchor_t::EvalBranch(std::vector<std::complex<double> >&branch, std::vector<std::complex<double> > &cpl, std::vector<double> &par, std::vector<double> &iso_par);

template<typename xdouble>
xdouble anchor_t::EvalAutoCpl(std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par,std::vector<xdouble> &iso_par){ // gets Chi2 with automatic non anchor couplings (no branchings)
	xdouble chi2 = 0.;
	int nNon = _borders_waves[0];
	std::vector<std::complex<xdouble> > actCpl = std::vector<std::complex<xdouble> >(nNon);
	for (int tbin=0;tbin<_nTbin;tbin++){
		if(_eval_tbin[tbin]){
			updateTprime(tbin);
			for(int i=0;i<nNon;i++){
				actCpl[i] = cpl[i+tbin*nNon];
			};
			chi2+=EvalAutoCplTbin(tbin,actCpl,par,iso_par);
	//		std::cout << tbin<<":A:"<< EvalAutoCplTbin(tbin,actCpl,par)<<std::endl;
		};
	};
	if (_write_out){
		*_outStream <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double anchor_t::EvalAutoCpl(std::vector<std::complex<double> > &cpl,std::vector<double> &par, std::vector<double> &iso_par);

template<typename xdouble>
xdouble anchor_t::EvalAutoCplBranch(std::vector<std::complex<xdouble> >&bra, std::vector<std::complex<xdouble> >&cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par){  // Gets Chi2 for automatically calculated couplings with branchings (the BEST)
	xdouble chi2=0.;
//	std::cout<<"par[1]: "<<par[1]<<std::endl;
	int nNon = _nBrCplAnc;
	std::vector<std::complex<xdouble> > actCpl = std::vector<std::complex<xdouble> >(nNon);
	for(int tbin=0;tbin<_nTbin;tbin++){
		if (_eval_tbin[tbin]){
			updateTprime(tbin);
			for(int i=0;i<nNon;i++){
				actCpl[i] = cpl[i+tbin*nNon];
			};
			std::vector<std::complex<xdouble> > bestcpl = getMinimumCplBra(tbin,bra,actCpl,par, iso_par);
		//	print_vector(bestcpl);
			std::vector<std::complex<xdouble> > best_cpl_std = std::vector<std::complex<xdouble> >(_nFtw);
			for (int i=0;i<_nFtw;i++){
				if(-1==_n_branch[i]){
					best_cpl_std[i] = bestcpl[_n_cpls[i]];
				}else{
					best_cpl_std[i] = bestcpl[_n_cpls[i]] * bra[_n_branch[i]]; 
				};
			};
			chi2+=EvalTbin(tbin,best_cpl_std,par,iso_par);
		};
	};
	if (_write_out){
		*_outStream <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double anchor_t::EvalAutoCplBranch(std::vector<std::complex<double> >&bra, std::vector<std::complex<double> >&cpl, std::vector<double> &par, std::vector<double> &iso_par);

template<typename xdouble>
xdouble anchor_t::EvalTbin(int tbin, std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par, std::vector<xdouble> &iso_par){ // gets the chi2 for a single t' bin
	updateTprime(tbin);
	xdouble chi2 = 0.;
	std::vector<std::vector<std::complex<xdouble> > > iso_eval = iso_funcs(iso_par);
	for (int bin=_minBin; bin<_maxBin; bin++){
		chi2+=EvalBin(tbin,bin,cpl,par,iso_eval);
//		std::cout<<"bin #"<<bin<<" chi2++"<<EvalBin(tbin,bin,cpl,par,iso_eval)<<std::endl;
//		std::cout << bin<<"   "<< EvalBin(tbin,bin,cpl,par) <<std::endl;
	};
	return chi2;
};
template double anchor_t::EvalTbin(int tbin, std::vector<std::complex<double> > &cpl,std::vector<double> &par, std::vector<double> &iso_par);

template<typename xdouble>
xdouble anchor_t::EvalAutoCplTbin(int tbin, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par){ // Gets the chi2 for a certain t' bin, with automatically calculated couplings
	std::vector<std::complex<xdouble> > mincpl = getMinimumCpl(tbin,cpl,par,iso_par);
	return EvalTbin(tbin,mincpl,par,iso_par);
};
template double anchor_t::EvalAutoCplTbin(int tbin, std::vector<std::complex<double> > &cpl, std::vector<double> &par, std::vector<double> &iso_par);

template<typename xdouble>
xdouble anchor_t::EvalBin(int tbin,int bin,std::vector<std::complex<xdouble> >&cpl, std::vector<xdouble> &par,std::vector<std::vector<std::complex<xdouble> > > &iso_eval){ // Gets the Chi2 for a single t' and m3pi bin
	double mass = (_binning[bin] + _binning[bin+1])/2; // Eval at bin center.
	std::vector<xdouble> deltas = delta(tbin,bin,mass, cpl,par,iso_eval);
	xdouble chi2 = 0.;
//	print_vector(deltas);
	for (int i=0;i<2*_nPoints-1;i++){
		int iWave = _point_to_wave[(i+1)/2];
		if (mass >= _lowerLims[iWave] and mass < _upperLims[iWave]){
//			std::cout<<"le_addite:D"<<pow(deltas[i],2.)*_coma[tbin][bin][i][i]<<std::endl;
			if (_write_out){
				*_outStream <<" mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<i+1<<"  ipj=            "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
				*_outStream <<"  delta1(imb, ipi, isectd) delta1(imb, ipj, isectd)     "<<deltas[i]<<"        "<<deltas[i]<<std::endl;
				*_outStream <<"   cov_key(imb, ipi, ipj,  isectd) =     "<<_coma[tbin][bin][i][i]<<std::endl;
				*_outStream <<"  Ci^2 before      "<<chi2<<"      +     "<<deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<"       =    ";
			};
			chi2+= deltas[i]*deltas[i]*_coma[tbin][bin][i][i];
//			std::cout<<deltas[i]<<"*"<<deltas[i]<<"*"<<_coma[tbin][bin][i][i]<<"="<<deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<std::endl;
			if(_write_out){
				*_outStream <<chi2<<std::endl;
			};
			for (int j=0;j<i;j++){
				int jWave = _point_to_wave[(j+1)/2];
				if(mass >= _lowerLims[jWave] and mass < _upperLims[jWave]){
//					std::cout<<"le_addite: "<<2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]<<std::endl;
					if (_write_out){
						*_outStream <<" mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<j+1<<"  ipj=            "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
						*_outStream <<"  delta1(imb, ipi, isectd) delta1(imb, ipj, isectd)     "<<deltas[j]<<"        "<<deltas[i]<<std::endl;
						*_outStream <<"   cov_key(imb, ipi, ipj,  isectd) =     "<<_coma[tbin][bin][j][i]<<std::endl;
						*_outStream <<"  Ci^2 before      "<<chi2<<"      +     "<<2.*deltas[j]*deltas[i]*_coma[tbin][bin][i][j]<<"       =    ";
					};
					chi2+=2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]; // Factor 2. because _coma[t][m][i][j] is symmetric under i<->j
//					std::cout<<deltas[i]<<"*"<<deltas[j]<<"*"<<_coma[tbin][bin][i][j]<<"="<<deltas[i]*deltas[j]*_coma[tbin][bin][i][j]<<std::endl;
					if(_write_out){
						*_outStream <<chi2<<std::endl;
					};
				};	
			};
		};
	};
	return chi2;
};
template double anchor_t::EvalBin(int tbin,int bin,std::vector<std::complex<double> >&cpl,std::vector<double> &par,std::vector<std::vector<std::complex<double> > > &iso_eval);

template<typename xdouble>
std::vector<xdouble> anchor_t::delta(int tbin, int bin,double mass, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par, std::vector<std::vector<std::complex<xdouble> > > &iso_eval){ // Returns f(m,...) - data[...] for each SDM entry in the fit
	std::vector<std::complex<xdouble> > ampls = amps(mass, cpl, par, iso_eval);
	std::vector<xdouble> del = std::vector<xdouble> (2*_nPoints-1);
	std::complex<xdouble> divider = std::complex<xdouble>(1.,0.); // When amplitudes are fitted, this is set to |ampls[0]|
	if(_is_ampl){
		divider = std::complex<xdouble>(pow(std::norm(ampls[0]),.5),0.);
	};
	if(_write_out){
		*_outStream << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           1  isectd=           "<<tbin+1<<std::endl;
		*_outStream << " Re1 data    "<<_data[tbin][bin][0]<<"     - theory   "<<std::norm(ampls[0]/divider)<<"       =    "<<std::norm(ampls[0]/divider) - _data[tbin][bin][0]<<std::endl;
	};
	del[0]=std::norm(ampls[0]/divider) - _data[tbin][bin][0];
	for (int i = 1; i<_nPoints;i++){
		std::complex<xdouble> inter = ampls[0]*std::conj(ampls[i])/divider;
		del[2*i-1]=real(inter) - _data[tbin][bin][2*i-1]; // real part
		del[2*i]=imag(inter) - _data[tbin][bin][2*i];    // imag part
		if(_write_out){
			*_outStream << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<2*i<<"  isectd=           "<<tbin+1<<std::endl;
			*_outStream << " Re data    "<<_data[tbin][bin][2*i-1]<<"     - theory   "<<real(inter)<<"       =    "<<real(inter) - _data[tbin][bin][2*i-1]<<std::endl;
			*_outStream << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<2*i+1<<"  isectd=           "<<tbin+1<<std::endl;
			*_outStream << " Im data    "<<_data[tbin][bin][2*i]<<"     - theory   "<<imag(inter)<<"       =    "<<imag(inter) - _data[tbin][bin][2*i]<<std::endl;
		};
	};
	return del;
};
template std::vector<double> anchor_t::delta(int tbin, int bin,double mass, std::vector<std::complex<double> > &cpl, std::vector<double> &par, std::vector<std::vector<std::complex<double> > > &iso_eval);

std::vector<std::vector<double> > anchor_t::getPlots(int tbin, std::vector<std::complex<double> >&cpl, std::vector<double> &par){ ///// DOES NOT WORK !!!!!
	updateTprime(tbin);
	std::vector<std::vector<double> > plots = std::vector<std::vector<double> >(2*_nPoints -1); // Store everything for one data set in one vector, with the structure [fit, data, error, fit, data, error, fit, ... ]
	std::vector<std::vector<std::complex<double> > > iso_eval; // If isobars were used
	for(int bin=0;bin<_nBins;bin++){
		double mass = (_binning[bin]+_binning[bin+1])/2.;
		std::vector<std::complex<double> > ampls = amps(mass, cpl ,par,iso_eval);
		std::complex<double> divider = std::complex<double>(1.,0.); // When amplitudes are fitted, this is set to |ampls[0]|
		if(_is_ampl){
			divider = std::complex<double>(pow(std::norm(ampls[0]),.5),0.);
		};
		plots[0].push_back(std::norm(ampls[0]/divider));
		plots[0].push_back(_data[tbin][bin][0]);
		plots[0].push_back(_coma[tbin][bin][0][0]); // Use simple diagonal errors here.
		for (int wave=1;wave<_nPoints;wave++){
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
template adouble anchor_t::EvalCP(std::vector<std::complex<adouble> > &cpl,std::vector<adouble> &par, std::vector<adouble> &iso_par);
template adouble anchor_t::EvalBranch(std::vector<std::complex<adouble> >&branch, std::vector<std::complex<adouble> > &cpl, std::vector<adouble> &par, std::vector<adouble> &iso_par);
template adouble anchor_t::EvalTbin(int tbin, std::vector<std::complex<adouble> > &cpl,std::vector<adouble> &par, std::vector<adouble> &iso_par);
template adouble anchor_t::EvalBin(int tbin,int bin,std::vector<std::complex<adouble> >&cpl,std::vector<adouble> &par, std::vector<std::vector<std::complex<adouble> > > &iso_par);
template std::vector<adouble> anchor_t::delta(int tbin, int bin,double mass, std::vector<std::complex<adouble> > &cpl, std::vector<adouble> &par, std::vector<std::vector<std::complex<adouble> > > &iso_eval);

template adouble anchor_t::EvalAutoCpl(std::vector<std::complex<adouble> > &cpl,std::vector<adouble> &par, std::vector<adouble> &iso_par);
template adouble anchor_t::EvalAutoCplBranch(std::vector<std::complex<adouble> >&bra, std::vector<std::complex<adouble> >&cpl, std::vector<adouble> &par, std::vector<adouble> &iso_par);
template adouble anchor_t::EvalAutoCplTbin(int tbin, std::vector<std::complex<adouble> > &cpl, std::vector<adouble> &par, std::vector<adouble> &iso_par);

template AandB<adouble> anchor_t::get_AB(int tbin,std::vector<std::complex<adouble> > &anchor_cpl, std::vector<adouble> &par);
template AandB<adouble> anchor_t::get_AB_iso(int tbin,std::vector<std::complex<adouble> > &anchor_cpl, std::vector<adouble> &par, std::vector<adouble> &iso_par);
template std::vector<std::complex<adouble> > anchor_t::getMinimumCpl(int tbin,std::vector<std::complex<adouble> > &anchor_cpl, std::vector<adouble> &par, std::vector<adouble> &iso_par);
template std::vector<std::complex<adouble> > anchor_t::getMinimumCplBra(int tbin, std::vector<std::complex<adouble> > &branch, std::vector<std::complex<adouble> > &anchor_cpl, std::vector<adouble> &par, std::vector<adouble> &iso_par);
#endif//ADOL_ON


std::string anchor_t::className(){
	return "anchor_t";
};




bool anchor_t::set_data(int tbin, int bin,std::vector<double> data){ // Set data points manually
	if (_data.size() != _nTbin){
		_data = std::vector<std::vector<std::vector<double> > >(_nTbin);
	};
	if (_data[tbin].size() != _nBins){
		_data[tbin] = std::vector<std::vector<double> >(_nBins);
	};
	_data[tbin][bin] = data;
	if (data.size() == 2*_nPoints -1){
		return true;
	};
	return false;
};

bool anchor_t::set_coma(int tbin, int bin, std::vector<std::vector<double> > coma){
//	print_matrix(coma);
	if(_coma.size() != _nTbin){
		_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_nTbin);
	};
	if(_coma[tbin].size() != _nBins){
		_coma[tbin] = std::vector<std::vector<std::vector<double> > >(_nBins);
	};
	_coma[tbin][bin] = coma;
	if (coma.size() == 2*_nPoints-1){
		int nerr = 0;
		for (unsigned int i=0;i<coma.size();i++){
			if (coma[i].size() != coma.size()){
				std::cout<<"Warning: Set coma is not quadratic."<<std::endl;
				return false;
			};
		};	
		return true;
	};
	return false;
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
	unsigned int nDat = 2*_nPoints-1;
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

int anchor_t::getNanc(){ 
	return _nBrCplAnc * _nTbin;
};

int anchor_t::getNtotAnc(){ 
	return 2*getNanc() + getNpar() + 2*getNbra() + getNiso();
};

void anchor_t::loadComa(int tbin, const char* comaFile){ // Loads the (inverted) covariance matrix for a t' bin from 'comaFile'
	_coma[tbin]=std::vector<std::vector<std::vector<double> > >();
	std::fstream data(comaFile,std::ios_base::in);
	double val;
	unsigned int nDat = 2*_nPoints-1;
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


void anchor_t::nullify(){ // Sets all data to 0.
	for (int tbin = 0; tbin < _nTbin;tbin++){
		for (int bin =0; bin<_nBins;bin++){
			int nDat = 2*_nPoints-1;
			for (int i=0;i<nDat;i++){
				_data[tbin][bin][i]=0.;	
			};	
		};
	};
};

void anchor_t::update_n_cpls(){
	waveset::update_n_cpls();
	_nBrCplAnc = _n_cpls[_borders_waves[0]-1]+1;
};


/*
		┌───────────────────────────────────────┐
		│	Es folgt eine Indexschlacht	│
		└───────────────────────────────────────┘
*/
// Maybe this is obsolete, since the thing subborting the de-isobarred case should do the same and not be significantly slower. Nevertheless, keep it at the moment
template<typename xdouble>
AandB<xdouble> anchor_t::get_AB(int tbin,std::vector<std::complex<xdouble> > &anchor_cpl, std::vector<xdouble> &par){// Gets the coefficients A and B for Chi2 = x^T*A*x + B*x + C  
	int nCplAnc = _borders_waves[0]; // Number of couplings for the anchor wave
	int nNon = _nFtw - nCplAnc;
	AandB<xdouble> AB(2*nNon);
	std::vector<xdouble> lines = std::vector<xdouble>(2*_nWaves-2,0.);  
	for(int bin=_minBin;bin<_maxBin;bin++){
		double m = (_binning[bin]+_binning[bin+1])/2;
		std::vector<std::complex<xdouble> > func = funcs(m,par);

		std::vector<double> phase = phase_space(m);
		std::complex<xdouble> ampAnc = std::complex<xdouble>(0.,0.);
		for(int i=0;i<nCplAnc;i++){
			ampAnc+= anchor_cpl[i]*func[_funcs_to_waves[i]]*phase[0];
		};
		if(_is_ampl){ // Divide ampAnc by |ampAnc| to get the amplitude with phase 0. // Need to be tested
			ampAnc/=pow(std::norm(ampAnc),.5);
		};
		xdouble RR = ampAnc.real();
		xdouble II = ampAnc.imag();
		for (int i =0;i<2*_nWaves-2;i++){ 
			xdouble val = (RR*RR+II*II - _data[tbin][bin][0])*_coma[tbin][bin][i+1][0];
			for (int j=1;j<2*_nWaves-1;j++){ 
				val+=-_data[tbin][bin][j]*_coma[tbin][bin][i+1][j];
			};
			lines[i]=val;			
		};
		xdouble intAnc = RR*RR+II*II;
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

			std::complex<xdouble> AB1 = ampAnc * conj(func[f1]) * phase[wave1];

			xdouble r1 = func[f1].real()*phase[wave1];
			xdouble i1 = func[f1].imag()*phase[wave1];
			//// (RR + j II)*(R - j I)*(r - j i) = (RR R r - RR I i - II I r + II R i)+j(II R r + II I i - RR r I - RR R i)

			xdouble RR1 = RR*r1+II*i1; // = Coefficient of the real part of the coupling 1 in the real part of wave one
			xdouble RI1 =-RR*i1-II*r1; // = Coefficient of the imag part of the coupling 1 in the real part of wave one
			xdouble IR1 = II*r1-RR*i1; // = Coefficient of the real part of the coupling 1 in the imag part of wave one
			xdouble II1 = II*i1-RR*r1; // = Coefficient of the imag part of the coupling 1 in the imag part of wave one

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
								

				std::complex<xdouble> AB2 = ampAnc * conj(func[f2]) * phase[wave2];

				xdouble r2 = func[f2].real()*phase[wave2];
				xdouble i2 = func[f2].imag()*phase[wave2];
				//// (RR + j II)*(R - j I)*(r - j i) = (RR R r - RR I i - II I r + II R i)+j(II R r + II I i - RR r I - RR R i)
				xdouble RR2 = RR*r2+II*i2;
				xdouble RI2 =-RR*i2-II*r2;
				xdouble IR2 = II*r2-RR*i2;
				xdouble II2 = II*i2-RR*r2;
//				std::cout<<"_coma[tbin][bin][2*wave1-1][2*wave2-1]: "<<_coma[tbin][bin][2*wave1-1][2*wave2-1]<<std::endl;
//				std::cout<< "AB1.real(): "<<AB1.real()<<std::endl;
//				std::cout<< "AB2.real(): "<<AB2.real()<<std::endl;

				AB.A[2*ind1  ][2*ind2  ]+= AB1.real() * _coma[tbin][bin][2*wave1-1][2*wave2-1] * AB2.real();
				AB.A[2*ind1  ][2*ind2  ]+= AB1.real() * _coma[tbin][bin][2*wave1-1][2*wave2  ] * AB2.imag();
				AB.A[2*ind1  ][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*wave1  ][2*wave2-1] * AB2.real();
				AB.A[2*ind1  ][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*wave1  ][2*wave2  ] * AB2.imag();

				AB.A[2*ind1  ][2*ind2+1]+=-AB1.real() * _coma[tbin][bin][2*wave1-1][2*wave2  ] * AB2.real();
				AB.A[2*ind1  ][2*ind2+1]+= AB1.real() * _coma[tbin][bin][2*wave1-1][2*wave2-1] * AB2.imag();
				AB.A[2*ind1  ][2*ind2+1]+=-AB1.imag() * _coma[tbin][bin][2*wave1  ][2*wave2  ] * AB2.real();
				AB.A[2*ind1  ][2*ind2+1]+= AB1.imag() * _coma[tbin][bin][2*wave1  ][2*wave2-1] * AB2.imag();

				AB.A[2*ind1+1][2*ind2  ]+=-AB1.real() * _coma[tbin][bin][2*wave1  ][2*wave2-1] * AB2.real();
				AB.A[2*ind1+1][2*ind2  ]+=-AB1.real() * _coma[tbin][bin][2*wave1  ][2*wave2  ] * AB2.imag();
				AB.A[2*ind1+1][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*wave1-1][2*wave2-1] * AB2.real();
				AB.A[2*ind1+1][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*wave1-1][2*wave2  ] * AB2.imag();

				AB.A[2*ind1+1][2*ind2+1]+= AB1.real() * _coma[tbin][bin][2*wave1  ][2*wave2  ] * AB2.real();
				AB.A[2*ind1+1][2*ind2+1]+=-AB1.real() * _coma[tbin][bin][2*wave1  ][2*wave2-1] * AB2.imag();
				AB.A[2*ind1+1][2*ind2+1]+=-AB1.imag() * _coma[tbin][bin][2*wave1-1][2*wave2  ] * AB2.real();
				AB.A[2*ind1+1][2*ind2+1]+= AB1.imag() * _coma[tbin][bin][2*wave1-1][2*wave2-1] * AB2.imag();
			};
		};
	};
	return AB;
};
template AandB<double> anchor_t::get_AB(int tbin,std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par);

template<typename xdouble>
AandB<xdouble> anchor_t::get_AB_iso(int tbin, std::vector<std::complex<xdouble> > &anchor_cpl,std::vector<xdouble> &par, std::vector<xdouble> &iso_par){
	int nCplAnc = _borders_waves[0]; // Number of couplings for the anchor wave		
	std::vector<std::vector<std::complex<xdouble> > > iso_eval  =  iso_funcs(iso_par);					//// ><>< Maybe just give pointer to already evaluated functions to save time???
	int nNon = _nFtw - nCplAnc;
	AandB<xdouble> AB(2*nNon);
	std::vector<xdouble> lines = std::vector<xdouble>(2*_nPoints-2,0.); 							//// >>>>  _nPoints
	for(int bin=_minBin;bin<_maxBin;bin++){
		double m = (_binning[bin]+_binning[bin+1])/2;
		std::vector<std::complex<xdouble> > func = funcs(m,par);

		std::vector<double> phase = phase_space(m);
		std::complex<xdouble> ampAnc = std::complex<xdouble>(0.,0.);
		for(int i=0;i<nCplAnc;i++){
			ampAnc+= anchor_cpl[i]*func[_funcs_to_waves[i]]*phase[0];
		};
		if(_is_ampl){ // Divide ampAnc by |ampAnc| to get the amplitude with phase 0. // Need to be tested
			ampAnc/=pow(std::norm(ampAnc),.5);
		};
		xdouble RR = ampAnc.real();
		xdouble II = ampAnc.imag();
		for (int i =0;i<2*_nPoints-2;i++){ 										//// >>>>  _nPoints 
			xdouble val = (RR*RR+II*II - _data[tbin][bin][0])*_coma[tbin][bin][i+1][0];
			for (int j=1;j<2*_nPoints-1;j++){ 									//// >>>>  _nPoints 
				val+=-_data[tbin][bin][j]*_coma[tbin][bin][i+1][j];
			};
			lines[i]=val;			
		};
		xdouble intAnc = RR*RR+II*II;
		int wave1 = 1;	// Start with 1, leave anchor wave out								//// |||| For de-isobarred, this stays, because anchor wave must be isobarred
		int up1 = _borders_waves[1];
		for (unsigned int nf1 = nCplAnc;nf1<_nFtw;nf1++){
			if (nf1==up1){
				wave1+=1;
				up1 = _borders_waves[wave1];
			};
			int f1 = _funcs_to_waves[nf1];
			int iso_f1 = _iso_to_waves[nf1];									//// >>>> Int f isobar 1
			int iso_nBin1 = abs(_wave_binning_pts[wave1]);
			int ind1=nf1-nCplAnc;
			int wave_start1 = _point_borders_wave[wave1-1];								//// >>>> To calculate the _data and _coma indices

			xdouble r1;
			xdouble i1;


			for (int i_iso_1=0;i_iso_1<iso_nBin1;i_iso_1++){								//// >>>> Here loop ofer isobar bins, if necessary (i_iso_1)(_wave_binning_pts[wave1)
				std::complex<xdouble> AB1;
				if (iso_f1 == -1){
					AB1 = ampAnc * conj(func[f1]) * phase[wave1];
					r1 = func[f1].real()*phase[wave1];
					i1 = func[f1].imag()*phase[wave1];
				}else{
					AB1 = ampAnc * conj(func[f1]*iso_eval[iso_f1][i_iso_1]) * phase[wave1];				//// >>>> Here take isobar function evaluation
					r1 = (func[f1]*iso_eval[iso_f1][i_iso_1]).real()* phase[wave1];					//// >>>> Here take isobar function evaluation
					i1 = (func[f1]*iso_eval[iso_f1][i_iso_1]).imag()* phase[wave1];					//// >>>> Here take isobar function evaluation
				};

				int point1 = wave_start1 + i_iso_1;	// Position of data and coma points
				int wave2=1;
				int up2=_borders_waves[1];
				//// (RR + j II)*(R - j I)*(r - j i) = (RR R r - RR I i - II I r + II R i)+j(II R r + II I i - RR r I - RR R i)

				xdouble RR1 = RR*r1+II*i1; // = Coefficient of the real part of the coupling 1 in the real part of wave one
				xdouble RI1 =-RR*i1-II*r1; // = Coefficient of the imag part of the coupling 1 in the real part of wave one
				xdouble IR1 = II*r1-RR*i1; // = Coefficient of the real part of the coupling 1 in the imag part of wave one
				xdouble II1 = II*i1-RR*r1; // = Coefficient of the imag part of the coupling 1 in the imag part of wave one

				AB.B[2*ind1  ] += AB1.real() * lines[2*point1-2] * 2;							//// >>>> now is Point
				AB.B[2*ind1  ] += AB1.imag() * lines[2*point1-1] * 2;							//// >>>> now is Point
		
				AB.B[2*ind1+1] +=-AB1.real() * lines[2*point1-1] * 2;							//// >>>> now is Point
				AB.B[2*ind1+1] += AB1.imag() * lines[2*point1-2] * 2;							//// >>>> now is Point

				for (unsigned int nf2= nCplAnc;nf2<_nFtw;nf2++){
					if(nf2==up2){
						wave2+=1;
						up2 = _borders_waves[wave2];
					};
					int f2 = _funcs_to_waves[nf2];
					int iso_f2 = _iso_to_waves[nf2];								//// >>>> Int f isobar 2
					int iso_nBin2 = abs(_wave_binning_pts[wave2]);
					int wave_start2 = _point_borders_wave[wave2-1];							//// >>>> To calculate the _data and _coma indices
					int ind2=nf2-nCplAnc;										//// |||| ind still applies!!!
					for (int i_iso_2=0;i_iso_2<iso_nBin2;i_iso_2++){						//// >>>> loop over second isobar bins (i_iso_2) (_wave_binning_pts[wave2)

						int point2 = wave_start2+i_iso_2;
						std::complex<xdouble> AB2;
						xdouble r2;
						xdouble i2;
				if (iso_f2 == -1){
					AB2 = ampAnc * conj(func[f2]) * phase[wave2];
					r2 = func[f2].real()*phase[wave2];
					i2 = func[f2].imag()*phase[wave2];
				}else{
					AB2 = ampAnc * conj(func[f2]*iso_eval[iso_f2][i_iso_2]) * phase[wave2];				//// >>>> Here take isobar function evaluation
					r2 = (func[f2]*iso_eval[iso_f2][i_iso_2]).real()* phase[wave2];					//// >>>> Here take isobar function evaluation
					i2 = (func[f2]*iso_eval[iso_f2][i_iso_2]).imag()* phase[wave2];					//// >>>> Here take isobar function evaluation
				};
						//// (RR + j II)*(R - j I)*(r - j i) = (RR R r - RR I i - II I r + II R i)+j(II R r + II I i - RR r I - RR R i)
						xdouble RR2 = RR*r2+II*i2;
						xdouble RI2 =-RR*i2-II*r2;
						xdouble IR2 = II*r2-RR*i2;
						xdouble II2 = II*i2-RR*r2;
		//				std::cout<<"_coma[tbin][bin][2*wave1-1][2*wave2-1]: "<<_coma[tbin][bin][2*wave1-1][2*wave2-1]<<std::endl;
		//				std::cout<< "AB1.real(): "<<AB1.real()<<std::endl;
		//				std::cout<< "AB2.real(): "<<AB2.real()<<std::endl;

						AB.A[2*ind1  ][2*ind2  ]+= AB1.real() * _coma[tbin][bin][2*point1-1][2*point2-1] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2  ]+= AB1.real() * _coma[tbin][bin][2*point1-1][2*point2  ] * AB2.imag();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*point1  ][2*point2-1] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*point1  ][2*point2  ] * AB2.imag();	//// >>>> are points1/2 now

						AB.A[2*ind1  ][2*ind2+1]+=-AB1.real() * _coma[tbin][bin][2*point1-1][2*point2  ] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2+1]+= AB1.real() * _coma[tbin][bin][2*point1-1][2*point2-1] * AB2.imag();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2+1]+=-AB1.imag() * _coma[tbin][bin][2*point1  ][2*point2  ] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2+1]+= AB1.imag() * _coma[tbin][bin][2*point1  ][2*point2-1] * AB2.imag();	//// >>>> are points1/2 now

						AB.A[2*ind1+1][2*ind2  ]+=-AB1.real() * _coma[tbin][bin][2*point1  ][2*point2-1] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2  ]+=-AB1.real() * _coma[tbin][bin][2*point1  ][2*point2  ] * AB2.imag();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*point1-1][2*point2-1] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*point1-1][2*point2  ] * AB2.imag();	//// >>>> are points1/2 now

						AB.A[2*ind1+1][2*ind2+1]+= AB1.real() * _coma[tbin][bin][2*point1  ][2*point2  ] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2+1]+=-AB1.real() * _coma[tbin][bin][2*point1  ][2*point2-1] * AB2.imag();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2+1]+=-AB1.imag() * _coma[tbin][bin][2*point1-1][2*point2  ] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2+1]+= AB1.imag() * _coma[tbin][bin][2*point1-1][2*point2-1] * AB2.imag();	//// >>>> are points1/2 now

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
						std::cout << AB1.real() <<" "<< _coma[tbin][bin][2*point1-1][2*point2-1] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.real() <<" "<< _coma[tbin][bin][2*point1-1][2*point2  ] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1  ][2*point2-1] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1  ][2*point2  ] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now

						std::cout <<-AB1.real() <<" "<< _coma[tbin][bin][2*point1-1][2*point2  ] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.real() <<" "<< _coma[tbin][bin][2*point1-1][2*point2-1] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
						std::cout <<-AB1.imag() <<" "<< _coma[tbin][bin][2*point1  ][2*point2  ] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1  ][2*point2-1] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now

						std::cout <<-AB1.real() <<" "<< _coma[tbin][bin][2*point1  ][2*point2-1] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout <<-AB1.real() <<" "<< _coma[tbin][bin][2*point1  ][2*point2  ] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1-1][2*point2-1] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1-1][2*point2  ] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now

						std::cout << AB1.real() <<" "<< _coma[tbin][bin][2*point1  ][2*point2  ] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout <<-AB1.real() <<" "<< _coma[tbin][bin][2*point1  ][2*point2-1] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
						std::cout <<-AB1.imag() <<" "<< _coma[tbin][bin][2*point1-1][2*point2  ] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1-1][2*point2-1] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
*/
					};
				};
			};
		};
	};
	return AB;
};
template AandB<double> anchor_t::get_AB_iso(int tbin,std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par, std::vector<double> &iso_par);

template<typename xdouble>
std::vector<std::complex<xdouble> > anchor_t::getMinimumCpl(int tbin,std::vector<std::complex<xdouble> > &anchor_cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par){	// Gets the bes couplings without branchings
	int nCplAnc = _borders_waves[0]; // Number of couplings for the anchor wave
	int nNon = _nFtw - nCplAnc;
	std::vector<std::complex<xdouble> > cpl = std::vector<std::complex<xdouble> >(_nFtw);
	for (int i=0; i<nCplAnc; i++){
		cpl[i] = anchor_cpl[i];
	};
	AandB<xdouble>AB(0);
	if (_has_isobars){
		AB = get_AB_iso(tbin,anchor_cpl,par,iso_par);
	}else{
		AB = get_AB(tbin,anchor_cpl,par);
	};

	std::vector<xdouble> D = cholesky::cholesky_solve(AB.A,AB.B);

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
		cpl[nCplAnc+i]=std::complex<xdouble>(-D[2*i]/2.,-D[2*i+1]/2.); 
	};
	return cpl;
};
template std::vector<std::complex<double> > anchor_t::getMinimumCpl(int tbin,std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par, std::vector<double> &iso_par);

template<typename xdouble>
std::vector<std::complex<xdouble> > anchor_t::getMinimumCplBra(int tbin, std::vector<std::complex<xdouble> > &branch, std::vector<std::complex<xdouble> > &anchor_cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par){ // Calculated all non anchor couplings (no need for fitting)
	if (0==_nBranch){ // Do not do the complicated stuff, when no branchings are used
		return getMinimumCpl(tbin,anchor_cpl,par, iso_par);
	}else{
		int nCplAnc = _borders_waves[0];
		int nNon = _nFtw - nCplAnc;
		std::vector<std::complex<xdouble> > cplAncBr = std::vector<std::complex<xdouble> >(nCplAnc);
		std::vector<std::complex<xdouble> > ccppll(_nBrCpl);
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
		AandB<xdouble>AB(0); // Here probably the if can be deleted and everything be done with 'get_AB_iso(...)' since it produces the same result, nevertheless, keep it atm
		if (_has_isobars){
			AB = get_AB_iso(tbin,cplAncBr,par,iso_par);
		}else{
			AB = get_AB(tbin,cplAncBr,par);
		};
		AandB<xdouble> ABprime(2*(_nBrCpl - _nBrCplAnc));
		for (int i =0;i<nNon;i++){ // Reshuffle the coefficients 
			int i_tot = i+nCplAnc;
			int i_cpl_tot = _n_cpls[i_tot];
			int i_cpl = i_cpl_tot - _nBrCplAnc;
			int i_br  = _n_branch[i_tot];
	//		std::cout<<"────────────────────────────────────────────────"<<std::endl;
	//		std::cout<<i<<" "<<i_tot<<" "<<i_cpl_tot<<" "<<i_cpl<<" "<<i_br<<std::endl;		
			if (-1 == i_br){ // If the funcion is uncoupled, just take the old coefficients
				ABprime.B[2*i_cpl  ]+= AB.B[2*i  ];
				ABprime.B[2*i_cpl+1]+= AB.B[2*i+1];
			}else if(i_cpl>-1){ // If the function is not coupled to the anchor wave, add it to the corresponding coupling
				ABprime.B[2*i_cpl  ]+= AB.B[2*i  ]*branch[i_br].real()   +   AB.B[2*i+1]*branch[i_br].imag();
				ABprime.B[2*i_cpl+1]+=-AB.B[2*i  ]*branch[i_br].imag()   +   AB.B[2*i+1]*branch[i_br].real();
			};
			for (int j=0; j<nNon;j++){
				int j_tot = j+nCplAnc;
				int j_cpl_tot = _n_cpls[j_tot];
				int j_cpl = j_cpl_tot -_nBrCplAnc;
				int j_br  = _n_branch[j_tot];
				if ((-1 == i_br or i_cpl>-1)and(-1 == j_br or j_cpl>-1)){
					std::complex<xdouble> ifak(1.,0.);
					std::complex<xdouble> jfak(1.,0.);
					if (i_br >-1){
						ifak = branch[i_br];
					};
					if (j_br >-1){
						jfak = branch[j_br];
					};
					ABprime.A[2*i_cpl  ][2*j_cpl  ]+=	 AB.A[2*i  ][2*j  ]*ifak.real()*jfak.real() +
										 AB.A[2*i  ][2*j+1]*ifak.real()*jfak.imag() +
										 AB.A[2*i+1][2*j  ]*ifak.imag()*jfak.real() +
										 AB.A[2*i+1][2*j+1]*ifak.imag()*jfak.imag() ;
					ABprime.A[2*i_cpl  ][2*j_cpl+1]+=	-AB.A[2*i  ][2*j  ]*ifak.real()*jfak.imag() +
										 AB.A[2*i  ][2*j+1]*ifak.real()*jfak.real() +
										-AB.A[2*i+1][2*j  ]*ifak.imag()*jfak.imag() +
										 AB.A[2*i+1][2*j+1]*ifak.imag()*jfak.real() ;
					ABprime.A[2*i_cpl+1][2*j_cpl  ]+=	-AB.A[2*i  ][2*j  ]*ifak.imag()*jfak.real() +
										-AB.A[2*i  ][2*j+1]*ifak.imag()*jfak.imag() +
										 AB.A[2*i+1][2*j  ]*ifak.real()*jfak.real() +
										 AB.A[2*i+1][2*j+1]*ifak.real()*jfak.imag() ;
					ABprime.A[2*i_cpl+1][2*j_cpl+1]+=	 AB.A[2*i  ][2*j  ]*ifak.imag()*jfak.imag() +
										-AB.A[2*i  ][2*j+1]*ifak.imag()*jfak.real() +
										-AB.A[2*i+1][2*j  ]*ifak.real()*jfak.imag() +
										 AB.A[2*i+1][2*j+1]*ifak.real()*jfak.real() ;
	//				std::cout<<ifak<<jfak<<std::endl;
				}else if(-1 == i_br or i_cpl>-1){			
					std::complex<xdouble> ifak(1.,0.);
					if (i_br>-1){
						ifak = branch[i_br];
					};
					std::complex<xdouble> cpl_j = branch[j_br]*anchor_cpl[j_cpl_tot];
					ABprime.B[2*i_cpl  ]+=   AB.A[2*i  ][2*j  ]*ifak.real()*cpl_j.real() +
								 AB.A[2*i  ][2*j+1]*ifak.real()*cpl_j.imag() +
								 AB.A[2*i+1][2*j  ]*ifak.imag()*cpl_j.real() +
								 AB.A[2*i+1][2*j+1]*ifak.imag()*cpl_j.imag() ;
					ABprime.B[2*i_cpl+1]+=	-AB.A[2*i  ][2*j  ]*ifak.imag()*cpl_j.real() +
								-AB.A[2*i  ][2*j+1]*ifak.imag()*cpl_j.imag() +
								 AB.A[2*i+1][2*j  ]*ifak.real()*cpl_j.real() +
								 AB.A[2*i+1][2*j+1]*ifak.real()*cpl_j.imag() ;
	//				std::cout<<"Nirrr"<<std::endl;
				}else if(-1 == j_br or j_cpl>-1){	
					std::complex<xdouble> jfak(1.,0.);
					if (j_br>-1){
						jfak = branch[j_br];
					};
					std::complex<xdouble> cpl_i = branch[i_br]*anchor_cpl[i_cpl_tot];
					ABprime.B[2*j_cpl  ]+=	 AB.A[2*i  ][2*j  ]*cpl_i.real()*jfak.real() +
								 AB.A[2*i  ][2*j+1]*cpl_i.real()*jfak.imag() +
								 AB.A[2*i+1][2*j  ]*cpl_i.imag()*jfak.real() +
								 AB.A[2*i+1][2*j+1]*cpl_i.imag()*jfak.imag() ;
					ABprime.B[2*j_cpl+1]+=	-AB.A[2*i  ][2*j  ]*cpl_i.real()*jfak.imag() +
								 AB.A[2*i  ][2*j+1]*cpl_i.real()*jfak.real() +
								-AB.A[2*i+1][2*j  ]*cpl_i.imag()*jfak.imag() +
								 AB.A[2*i+1][2*j+1]*cpl_i.imag()*jfak.real() ;
	//				std::cout<<"Njrrr"<<std::endl;
				};// If both functions are coupled to the anchor wave, they just add to the constant C, (Chi2 = x^TAx + B^Tx +C) and do not change the position of the minimum.
			};
		};
	//	std::vector<double> les_as;
	//	for (int i=0;i<2*(_nBrCpl - _nBrCplAnc);i++){
	//		les_as.push_back(ABprime.A[i][i]);
	//	};
	//	print_vector(les_as);

		std::vector<xdouble> D = cholesky::cholesky_solve(ABprime.A,ABprime.B);

		for (int i=_nBrCplAnc;i<_nBrCpl;i++){
			int irel = i-_nBrCplAnc;
			ccppll[i] = std::complex<xdouble>(-D[2*irel]/2,-D[2*irel+1]/2);
		};
	//	std::cout<<"made it through"<<std::endl;
		return ccppll;
	};
};
template std::vector<std::complex<double> > anchor_t::getMinimumCplBra(int tbin, std::vector<std::complex<double> > &branch, std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par, std::vector<double> &iso_par);


void anchor_t::printStatus(){ // Prints the internal status
	waveset::printStatus();
	std::cout<<std::endl<< "_is_ampl: "<<_is_ampl<<std::endl;
	std::cout<<std::endl<< "_nBins: "<<_nBins<<std::endl;
	std::cout<<std::endl<< "_minBin: "<<_minBin<<std::endl;
	std::cout<<std::endl<< "_maxBin: "<<_maxBin<<std::endl;
	std::cout<<std::endl<< "_nTbin: "<<_nTbin<<std::endl;
	std::cout<<std::endl<< "_mMin: "<<_mMin<<std::endl;
	std::cout<<std::endl<< "_mMax: "<<_mMax<<std::endl;
	std::cout<<std::endl<<"_const_is_t:"<<std::endl;
	print_vector(_const_is_t);
	std::cout<<std::endl<<"_t_binning:"<<std::endl;
	print_vector(_t_binning);
	std::cout<<std::endl<<"_binning:"<<std::endl;
	print_vector(_binning);
	std::cout<<std::endl<<"Omit printing of _data and _coma"<<std::endl;
	std::cout<<std::endl<< "_nBranch: "<<_nBranch<<std::endl;
	std::cout<<std::endl<< "_nBrCpl: "<<_nBrCpl<<std::endl;
	std::cout<<std::endl<< "_nBrCplAnc: "<<_nBrCplAnc<<std::endl;
	std::cout<<std::endl<< "_coupled:"<<std::endl;
	print_vector(_coupled);
	std::cout<<std::endl<< "_n_branch:"<<std::endl;
	print_vector(_n_branch);
	std::cout<<std::endl<< "_n_cpls:"<<std::endl;	
	print_vector(_n_cpls);
	std::cout<<std::endl<< "_eval_tbin"<<std::endl;
	print_vector(_eval_tbin);	
};




void anchor_t::set_is_ampl(bool is_ampl){ // Set the mode to amplitude fitting, instead of anchor wave
	_is_ampl = is_ampl;
};

std::vector<std::complex<double> > anchor_t::get_branchings(std::vector<std::complex<double> > &cpl,std::vector<double> &par, std::vector<double> &iso_par){ 
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
			std::vector<std::complex<double> > all_cpls_t= getMinimumCpl(tbin,cpl_t,par, iso_par);
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

////////////////////////////// WILL STAY!!!!!!!!!!!!!!!!!!!!!!


void anchor_t::conjugate(){
	for (int i=2;i<2*(_nPoints);i+=2){
		for (int tbin=0;tbin<_nTbin;tbin++){
			for(int bin=0;bin<_nBins;bin++){
				_data[tbin][bin][i]*=-1;
				for (int j=0;j<2*_nPoints-1;j++){
					_coma[tbin][bin][i][j]*=-1;
					_coma[tbin][bin][j][i]*=-1;
				};
			};
		};
	};
};

void anchor_t::setTbinning(std::vector<double> binning){ /// Has to stay, gives right size for _data and _coma
	waveset::setTbinning(binning);
	if (_data.size() != _nTbin){
		if (_data.size() != 0){
			std::cout<<"Warning: _data.size() != _nTbin, but nonzero. Previous set _data[][][] will be lost."<<std::endl;
		};
		_data = std::vector<std::vector<std::vector<double> > >(_nTbin);
	};
	if (_coma.size() != _nTbin){
		if (_coma.size() != 0){
			std::cout<<"Warning: _coma.size() != _nTbin, bun nonzero. Previous set _coma[][][][] will be lost."<<std::endl;
		};
		_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_nTbin);
	};
};
