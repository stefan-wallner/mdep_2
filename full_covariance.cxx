#include"full_covariance.h"
#include<vector>
#include<complex>
#include<fstream>
#include<iostream>
#include<string>
#include<limits>
#include <stdexcept>
#include"matrix_utilities.h"

#ifdef ADOL_ON
#include "adolc/adolc.h"
std::complex<adouble> operator*(std::complex<adouble> a,double d){// Define std::complex<adouble> * double
	return std::complex<adouble>(a.real()*d,a.imag()*d);
};
#endif//ADOL_ON

full_covariance::full_covariance():
	_waveset(),
	_nOut(1000),
	_count(0),
	_nBra(0){
	
	setTbinning((*_waveset.t_binning()));
};
//########################################################################################################################################################
#ifdef USE_YAML
///Constructror from YAML files
full_covariance::full_covariance(
							std::string 						card):

	_waveset(card),
	_nOut(1000),
	_count(0),
	_nBra(0){
	YAML::Node Ycard   = YAML::LoadFile(card);
	std::string parametrizations = Ycard["parametrization_file"].as<std::string>();
	YAML::Node Yparam  = YAML::LoadFile(parametrizations);
	update_n_cpls();
	setTbinning((*_waveset.t_binning()));
	update_definitions();
	update_min_max_bin();
	std::cout<<"Load full_covariance from YAML file\nLoad data and coma"<<std::endl;
	loadDataComa(Ycard);
	std::cout<<"Data and coma loaded\nLoad parameter values"<<std::endl;
	loadParameterValues(Ycard, Yparam);
	std::cout<<"Paramter values loaded"<<std::endl;
	std::cout<<"full_covariance loaded"<<std::endl;
};
#endif//USE_YAML
//########################################################################################################################################################
///Evaluate Chi2 with the internal _parameters (Call the other operator)
double full_covariance::operator()(){

	return (*this)(&parameters()[0]);
};
//########################################################################################################################################################
///Evaluate Chi2 with the paramters xx
double full_covariance::operator()(
							std::vector<double> 						&xx){

	return (*this)(&xx[0]);
};
//########################################################################################################################################################
///Evaluate Chi2 with the paramters xx
double full_covariance::operator()(
							const double							*xx){


//	const std::complex<double>* cpl = (std::complex<double>*)xx; // This is forbidden by Charly, build complex variables by hand :(
	std::complex<double> cpl[_nCpl];
	for (size_t i=0;i<_nCpl;i++){
		cpl[i] = std::complex<double>(xx[2*i],xx[2*i+1]);
	};
	const double* par = xx+2*_nCpl;
//	const std::complex<double>* bra = (std::complex<double>*)(xx+2*_nCpl+_nPar);// This is forbidden by Charly, build complex variables by hand :(
	std::complex<double> bra[_nBra];
	for (size_t i=0;i<_nBra;i++){
		bra[i] = std::complex<double>(xx[2*_nCpl+_nPar+2*i],xx[2*_nCpl+_nPar+2*i+1]);
	};
	const double* iso_par = xx + 2*_nCpl + _nPar + 2*_nBra;

	// std::cout<<par[0]<<std::endl;
	double chi2;
	chi2 = EvalBranch(bra,cpl,par,iso_par);
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
//########################################################################################################################################################
///Evaluates chi2 w/o branchings (EvalC(ouplings)P(arameters))
template<typename xdouble>
xdouble full_covariance::EvalCP(
							const std::complex<xdouble>					*cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)		const{

	xdouble chi2 = 0.;
	std::vector<std::complex<xdouble> > actCpl = std::vector<std::complex<xdouble> >(_waveset.nFtw());
	for (size_t tbin=0;tbin<_waveset.nTbin();tbin++){
		if((*_waveset.eval_tbin())[tbin]){
			for (size_t i=0;i<_waveset.nFtw();i++){
				actCpl[i] = cpl[i+tbin*_waveset.nFtw()];
			};
//			std::cout<<"VVVVVVVVVVVVVVVVVVVVVVVVV"<<std::endl;
//			print_vector(actCpl);
//			std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
			chi2+=EvalTbin(tbin,&actCpl[0],&par[0],&iso_par[0]);//[0]//
	//		std::cout<<tbin<<":E:"<<EvalTbin(tbin,actCpl,par)<<std::endl;
//			std::cout<<chi2<<std::endl;
		};
	};
	if (_waveset.write_out()){
		*_waveset.outStream() <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double full_covariance::EvalCP(const std::complex<double> *cpl,const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Evaluates chi2 with branchings
template<typename xdouble>
xdouble full_covariance::EvalBranch(
							const std::complex<xdouble>					*branch,
							const std::complex<xdouble>	 				*cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)		const{

	std::vector<std::complex<xdouble> > cplin(_waveset.nFtw()*_waveset.nTbin());
	for (size_t tbin=0;tbin<_waveset.nTbin();tbin++){
		if((*_waveset.eval_tbin())[tbin]){
			for (size_t i =0;i<_waveset.nFtw();i++){
				if ((*_waveset.n_branch())[i]==-1){
					cplin[i+tbin*_waveset.nFtw()] = cpl[(*_waveset.n_cpls())[i]+tbin*(_waveset.nBrCpl())];
				}else{
					cplin[i+tbin*_waveset.nFtw()] = cpl[(*_waveset.n_cpls())[i]+tbin*_waveset.nBrCpl()] * branch[(*_waveset.n_branch())[i]];
				};
			};
		};
	};
	xdouble chi2 = EvalCP(&cplin[0],par, iso_par);//&...[0] stays
	if (_waveset.write_out()){
		*_waveset.outStream() <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double full_covariance::EvalBranch(const std::complex<double> *branch, const std::complex<double> *cpl, const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Gets the chi2 for a single t' bin
double full_covariance::EvalTbin(
							int 								tbin,
							const std::complex<double>	 				*cpl,
							const double	 						*par,
							const double 							*iso_par)		const{

	double chi2 = 0.;
	std::vector<std::vector<std::complex<double> > > iso_eval = _waveset.iso_funcs(iso_par);
//Careful, multiprocessing messes up the output file!!! Bed for debugging
#pragma omp parallel for reduction(+:chi2)
	for (size_t bin=_waveset.minBin(); bin<_waveset.maxBin(); bin++){
		chi2+=EvalBin(tbin,bin,cpl,par,iso_eval);
//		std::cout<<"bin #"<<bin<<" chi2++"<<EvalBin(tbin,bin,cpl,par,iso_eval)<<std::endl;
//		std::cout << bin<<"   "<< EvalBin(tbin,bin,cpl,par) <<std::endl;
	};
	return chi2;
};
#ifdef ADOL_ON // Has to be twice, otherwise the openmp pragma would not work. Did not fins another workaround
adouble full_covariance::EvalTbin(
							int 								tbin,
							const std::complex<adouble>	 				*cpl,
							const adouble	 						*par,
							const adouble 							*iso_par)		const{

	adouble chi2 = 0.;
	std::vector<std::vector<std::complex<adouble> > > iso_eval = _waveset.iso_funcs(iso_par);
	for (size_t bin=_waveset.minBin(); bin<_waveset.maxBin(); bin++){
		chi2+=EvalBin(tbin,bin,cpl,par,iso_eval);
//		std::cout<<"bin #"<<bin<<" chi2++"<<EvalBin(tbin,bin,cpl,par,iso_eval)<<std::endl;
//		std::cout << bin<<"   "<< EvalBin(tbin,bin,cpl,par) <<std::endl;
	};
	return chi2;
};
#endif//ADOL_ON
//########################################################################################################################################################
///Gets the Chi2 for a single t' and m3pi bin
template<typename xdouble>
xdouble full_covariance::EvalBin(
							int 								tbin,
							int 								bin,
							const std::complex<xdouble>					*cpl,
							const xdouble	 						*par,
							std::vector<std::vector<std::complex<xdouble> > > 		&iso_eval)		const{

	double mass = _waveset.get_m(bin); // Eval at bin center.

//	std::vector<std::complex<xdouble> > ccpl(_waveset.nFtw(),std::complex<xdouble>(1.,0.));
//	std::cout<<ccpl.size()<<std::endl;
//	std::cout<<_nCpl<<std::endl;
//	for (size_t i=0; i<_nPar;++i){
//		std::cout<<i<<" "<<par[i]<<std::endl;
//	};


	std::vector<xdouble> deltas = delta(tbin,bin,mass, cpl, par,iso_eval);


	xdouble chi2 = 0.;
//	print_vector(deltas);
	for (size_t i=0;i<_waveset.nPoints()*_waveset.nPoints();++i){
//		int iWave1 = (*_waveset.point_to_wave())[(i+1)/2];
//		if (mass >= _lowerLims[iWave] and mass < _upperLims[iWave]){
//		if (_is_active[bin][i]){
		if (deltas[i] != 0.){

//			std::cout<<"le_addite:D"<<pow(deltas[i],2.)*_coma[tbin][bin][i][i]<<std::endl;
			if (_waveset.write_out()){
				*_waveset.outStream() <<" mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<i+1<<"  ipj=            "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
				*_waveset.outStream() <<"  delta1(imb, ipi, isectd) delta1(imb, ipj, isectd)     "<<deltas[i]<<"        "<<deltas[i]<<std::endl;
				*_waveset.outStream() <<"   cov_key(imb, ipi, ipj,  isectd) =     "<<_coma[tbin][bin][i][i]<<std::endl;
				*_waveset.outStream() <<"  Ci^2 before      "<<chi2<<"      +     "<<deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<"       =    ";
			};

//			if (tbin == 3 and bin == 50){
//				std::cout<< i<<" "<<deltas[i]<<" "<<_coma[tbin][bin][i][i]<<std::endl;
//			};

			chi2+= deltas[i]*deltas[i]*_coma[tbin][bin][i][i];
//			std::cout<<deltas[i]<<"*"<<deltas[i]<<"*"<<_coma[tbin][bin][i][i]<<"="<<deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<std::endl;
			if(_waveset.write_out()){
				*_waveset.outStream() <<chi2<<std::endl;
			};

			for (size_t j=0;j<i;j++){
//				int jWave = (*_waveset.point_to_wave())[(j+1)/2];
//				if(mass >= _lowerLims[jWave] and mass < _upperLims[jWave]){
//				if(_is_active[bin][j]){
				if (deltas[j] != 0.){
//					std::cout<<"le_addite: "<<2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]<<std::endl;
					if (_waveset.write_out()){
						*_waveset.outStream() <<" mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<j+1<<"  ipj=            "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
						*_waveset.outStream() <<"  delta1(imb, ipi, isectd) delta1(imb, ipj, isectd)     "<<deltas[j]<<"        "<<deltas[i]<<std::endl;
						*_waveset.outStream() <<"   cov_key(imb, ipi, ipj,  isectd) =     "<<_coma[tbin][bin][j][i]<<std::endl;
						*_waveset.outStream() <<"  Ci^2 before      "<<chi2<<"      +     "<<2.*deltas[j]*deltas[i]*_coma[tbin][bin][i][j]<<"       =    ";
					};
					chi2+=2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]; // Factor 2. because _coma[t][m][i][j] is symmetric under i<->j
//					std::cout<<deltas[i]<<"*"<<deltas[j]<<"*"<<_coma[tbin][bin][i][j]<<"="<<deltas[i]*deltas[j]*_coma[tbin][bin][i][j]<<std::endl;
					if(_waveset.write_out()){
						*_waveset.outStream() <<chi2<<std::endl;
					};
				};
			};
		};
	};
//	std::cout<<tbin<<":"<<bin<<"::"<<chi2<<std::endl;
	return chi2;
};
template double full_covariance::EvalBin(int tbin,int bin,const std::complex<double> *cpl,const double *par,std::vector<std::vector<std::complex<double> > > &iso_eval) const;
//########################################################################################################################################################
///Returns f(m,...) - data[...] for each SDM entry in the fit
template<typename xdouble>
std::vector<xdouble> full_covariance::delta(
							int 								tbin,
							int 								bin,
							double 								mass,
							const std::complex<xdouble> 					*cpl,
							const xdouble	 						*par,
							std::vector<std::vector<std::complex<xdouble> > > 		&iso_eval)		const{

	size_t nPoints = _waveset.nPoints();
	std::vector<double> var = _waveset.getVar(mass,tbin);
	std::vector<std::complex<xdouble> > ampls = _waveset.amps(&var[0], cpl, par, iso_eval);


///	std::cout<<"''''''''''''''''''''''''''''''''''''''''"<<std::endl;
///	for (size_t i=0;i<ampls.size();++i){
///		std::cout<<i<<": "<<ampls[i]<<std::endl;
///	};
	std::vector<xdouble> del(nPoints*nPoints,0.);
//	del[0]=std::norm(ampls[0]) - _data[tbin][bin][0]; ////// HERE ADD if(_is_point_bin[bin][i]){...}; so that no deltas are aclculated for turend off bins
	for (size_t i = 0; i<_waveset.nPoints();i++){
//		if (not _is_active[bin][i]){

//			continue;
//		};
		if (ampls[i] == std::complex<xdouble>(0.,0.)){
			continue;
		};

		del[i*(nPoints+1)] = std::norm(ampls[i]) - _data[tbin][bin][i*(nPoints+1)];
		if(_waveset.write_out()){
			*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
			*_waveset.outStream() << " Re1 data    "<<_data[tbin][bin][i*(nPoints+1)]<<"     - theory   "<<std::norm(ampls[i])<<"       =    "<<del[i*(nPoints+1)]<<std::endl;
		};
		for (size_t j=0; j<i;++j){
//			if (not _is_active[bin][j]){
//				continue;
//			};
			
			if (ampls[j] == std::complex<xdouble>(0.,0.)){
				continue;
			};

			std::complex<xdouble> inter = ampls[i]*std::conj(ampls[j]);
			del[nPoints*i+j]=   imag(inter) - _data[tbin][bin][nPoints*i+j];    // imag part
			del[nPoints*j+i]=   real(inter) - _data[tbin][bin][nPoints*j+i];    // real part

			if(_waveset.write_out()){
				*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<nPoints*i+j<<"  isectd=           "<<tbin+1<<std::endl;
				*_waveset.outStream() << " Re data    "<<_data[tbin][bin][nPoints*i+j]<<"     - theory   "<<real(inter)<<"       =    "<<del[nPoints*i+j]<<std::endl;
				*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<nPoints*j+i<<"  isectd=           "<<tbin+1<<std::endl;
				*_waveset.outStream() << " Im data    "<<_data[tbin][bin][nPoints*j+i]<<"     - theory   "<<imag(inter)<<"       =    "<<del[nPoints*j+i]<<std::endl;
			};
		};
	};
//	print_vector(del);	

	return del;
};
template std::vector<double> full_covariance::delta(int tbin, int bin,double mass, const std::complex<double> *cpl, const double *par, std::vector<std::vector<std::complex<double> > > &iso_eval) const;
//########################################################################################################################################################
///Instantiate auto diff methods, if needed (Enable adouble operations, if the auto diff package is loaded)
#ifdef ADOL_ON
template adouble full_covariance::EvalCP(const std::complex<adouble> *cpl,const adouble *par, const adouble *iso_par) const;
template adouble full_covariance::EvalBranch(const std::complex<adouble> *branch, const std::complex<adouble> *cpl, const adouble *par, const adouble *iso_par) const;
template adouble full_covariance::EvalBin(int tbin,int bin,const std::complex<adouble> *cpl,const adouble *par, std::vector<std::vector<std::complex<adouble> > > &iso_par) const;
template std::vector<adouble> full_covariance::delta(int tbin, int bin,double mass, const std::complex<adouble> *cpl, const adouble *par, std::vector<std::vector<std::complex<adouble> > > &iso_eval) const;
//#######################################################################################################################################################
///Gets the gradient w.r.t. xx
std::vector<double> full_covariance::Diff(
							std::vector<double> 				&xx)						const{

	return Diff(&xx[0]);
};
//#######################################################################################################################################################
///Gets the gradient w.r.t. xx
std::vector<double> full_covariance::Diff(
							const double 					*xx)						const{

	int nTape = 0;
	double x[_nTot];
	for (size_t i=0;i<_nTot;i++){
		x[i] = xx[i];
	};

	trace_on(nTape);
	std::vector<adouble>aCpl_r(2*_nCpl);
	std::vector<adouble>aPar(_nPar);
	std::vector<adouble>aBra_r(2*_nBra);
	std::vector<adouble>aIso(_nIso);
	int count=0;
	for (size_t i=0;i<2*_nCpl;i++){
		aCpl_r[i] <<= x[count];
		++count;
	};
	for (size_t i=0; i<_nPar;i++){
		aPar[i] <<= x[count];
		++count;
	};
	for (size_t i=0;i<2*_nBra;i++){
		aBra_r[i] <<= x[count];
		++count;
	};
	for (size_t i=0;i<_nIso;i++){
		aIso[i] <<= x[count];
		++count;
	};
	std::vector<std::complex<adouble> > aCpl_c(_nCpl);
	std::vector<std::complex<adouble> > aBra_c(_nBra);
	for (size_t i=0;i<_nCpl;i++){
		aCpl_c[i] = std::complex<adouble>(aCpl_r[2*i],aCpl_r[2*i+1]);
	};
	for (size_t i=0;i<_nBra;i++){
		aBra_c[i] = std::complex<adouble>(aBra_r[2*i],aBra_r[2*i+1]);
	};
	double Chi2;
	adouble aChi2;
	aChi2 = EvalBranch(&aBra_c[0],&aCpl_c[0],&aPar[0],&aIso[0]);//[0]//
	aChi2 >>= Chi2;
	trace_off();
	double grad[_nTot];
	gradient(nTape,_nTot,x,grad);
	vector<double> gradient(_nTot);
	for (size_t i=0;i<_nTot;i++){
		gradient[i]=grad[i];
	};
	return gradient;
};
#endif//ADOL_ON
//#######################################################################################################################################################
///Set a single parameter by number
void full_covariance::setParameter(
							size_t 						i, 	// # of parameter
							double 						par){

	if(_parameters.size() < 2*(_nCpl+_nBra)){
		_parameters = std::vector<double>(2*(_nCpl+_nBra),0.);
	};
	if (i<2*_nCpl){
		_parameters[i] = par;
	}else if(i<2*_nCpl+_nPar){
		_waveset.setPar(i-2*_nCpl,par);
	}else if(i<2*_nCpl+_nPar+2*_nBra){
		_parameters[i-_nPar] = par;
	}else if(i<2*_nCpl+_nPar+2*_nBra+_nIso){
		_waveset.setIsoPar(i-2*_nCpl-_nPar-2*_nBra,par);
	}else{
		std::cerr<<"Error: Can't set parameter #"<<i<<std::endl;
	};
};
//#######################################################################################################################################################
double full_covariance::getParameter(			size_t 						i)						const{

	if(i<2*_nCpl){
		return _parameters[i];
	}else if(i<2*_nCpl+_nPar){
		return _waveset.getPar(i-2*_nCpl);
	}else if(i<2*_nCpl+_nPar+2*_nBra){
		return _parameters[i-_nPar];
	}else if(i<2*_nCpl+_nPar+2*_nBra+_nIso){
		return _waveset.getIsoPar(i-2*_nCpl-_nPar-2*_nBra);
	}else{
		std::cerr<<"Error: Can't get parameter #"<<i<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};
//#######################################################################################################################################################
///Set internal parameters
void full_covariance::setParameters(
							std::vector<double> 				pars){

	if(pars.size()>=_nTot){
		for (size_t i=0;i<_nTot;i++){
			setParameter(i,pars[i]);
		};
	}else{
		std::cerr<<"Error: Input pars.size() too small"<<std::endl;
	};
};
//#######################################################################################################################################################
///Get parameters
const std::vector<double> full_covariance::parameters()													const{

	std::vector<double> par(_nTot);
	for (size_t i=0;i<_nTot;i++){
		par[i] = getParameter(i);
	};
	return par;
};
//#######################################################################################################################################################
///Set a single parameter by name
bool full_covariance::setParameter(
							std::string 					name,
							double 						par){

	int number = getParNumber(name);
	if (-1==number){
		std::cerr << "Error: Parameter '"<<name<<"' not found"<<std::endl;
		return false;
	}else if ((int)_nTot <= number){
		std::cerr << "Error: Parameter number too high"<<std::endl;
		return false;
	}else{
		setParameter(number,par);
		return true;
	};
};
//########################################################################################################################################################
///Gets the number for a given name, if the name does not exist, return -1
int full_covariance::getParNumber(
							std::string 					name)						const{

	for (size_t i=0;i<_parNames.size();i++){
		if (_parNames[i] == name){
			return i;
		};
	};
	return -1;
};
//########################################################################################################################################################
///Set Parameter limits
void full_covariance::setParLimits(
							int 						par,
							double 						upper,
							double 						lower){

	if(_lower_parameter_limits.size() != _nTot){
		init_lower_limits();
	};
	if(_upper_parameter_limits.size() != _nTot){
		init_upper_limits();
	};
	_upper_parameter_limits[par]=upper;
	_lower_parameter_limits[par]=lower;
};
//########################################################################################################################################################
///Inits all lower parateter limits to 1.
void full_covariance::init_lower_limits(int n){

	if (n==-1){
		n=_nTot;
	};
	_lower_parameter_limits = std::vector<double>(n,1.);
};
//########################################################################################################################################################
///Inits all upper parateter limits to 0.
void full_covariance::init_upper_limits(int n){

	if(n==-1){
		n=_nTot;
	};
	_upper_parameter_limits = std::vector<double>(n,0.);
};
//########################################################################################################################################################
///Calculates some initial values for the branchigs
std::vector<std::complex<double> > full_covariance::get_branchings(
							const std::vector<std::complex<double> > 	&cpl,
							const std::vector<double> 			&par = std::vector<double>(),
							const std::vector<double> 			&iso_par = std::vector<double>())		const{

	


	throw; // Does not work at the moment

	return std::vector<std::complex<double> >(_nBra,std::complex<double>(1.,0.));;
};
//########################################################################################################################################################
///Gets the couplings-without-branchings from couplings-with-branchings and branchings for one t' bin
std::vector<std::complex<double> > full_covariance::getUnbranchedCouplings(
							const std::vector<std::complex<double> >	&cpl,
							const std::vector<std::complex<double> >	&bra)						const{
	
	std::vector<std::complex<double> > cpl_un(_waveset.nFtw());
	for (size_t i=0;i<_waveset.nFtw();i++){
		if(-1==(*_waveset.n_branch())[i]){
			cpl_un[i] = cpl[(*_waveset.n_cpls())[i]];
		}else{
			cpl_un[i] = cpl[(*_waveset.n_cpls())[i]] * bra[(*_waveset.n_branch())[i]];
		};
	};
	return cpl_un;
};
//########################################################################################################################################################
///Gets all couplings for one t' bin (.size() = _waveset.nFtw()) for different input formats
std::vector<std::complex<double> > full_covariance::getAllCouplings(
							int						tbin,
							const std::vector<std::complex<double> >	&cpl,
							const std::vector<double>			&par,
							const std::vector<std::complex<double> >	&bra,
							const std::vector<double>			&iso)						const{


	size_t nTbin = _waveset.nTbin();
	size_t nFtw = _waveset.nFtw();
	size_t cpls = _nCpl/nTbin;
	std::vector<std::complex<double> > cpl_all(nFtw);
	if (cpl.size() == nFtw*nTbin){
		std::cout<<"Got all couplings for all 't bins"<<std::endl;
		for (size_t i=0;i<nFtw;++i){
			cpl_all[i] = std::complex<double>(cpl[nFtw*tbin+i]);
		};	
	} else if (cpl.size() == _nCpl){
		std::cout<< "Got branched couplings for all t' bins"<<std::endl;
		std::vector<std::complex<double> > tCpls(cpls);
		for (size_t i=0;i<cpls;++i){
			tCpls[i]=cpl[tbin*cpls+i];
		};
		cpl_all = getUnbranchedCouplings(tCpls,bra);
	} else if (cpl.size() == nFtw){
		std::cout<<"Got all couplings for on t' bin"<<std::endl;
		cpl_all = cpl;
	} else if (cpl.size() == cpls){
		std::cout<<"Got branched couplings for one t' bin"<<std::endl;
		cpl_all = getUnbranchedCouplings(cpl,bra);
	}else{
		throw std::invalid_argument("Unknown format for couplings");
	};


	return cpl_all;
};
//########################################################################################################################################################
///Set the anchor coupligs to one, that are coupled to another function
void full_covariance::branchCouplingsToOne(){
	for (size_t i=0;i<(*_waveset.coupled()).size();i++){
		if((*_waveset.coupled())[i] >=0){
			size_t nCpl = (*_waveset.n_cpls())[i];
			if (nCpl <= _nBrCpl){
				for (size_t tbin =0;tbin<_waveset.nTbin();tbin++){
					setParameter(2*_nBrCpl*tbin+2*nCpl  ,1.);
					setParameter(2*_nBrCpl*tbin+2*nCpl+1,0.);
				};
			};
		};
	};
};
//########################################################################################################################################################
///Set data points manually
bool full_covariance::set_data(
							int								tbin,
							int 								bin,
							std::vector<double> 						data){

	// Ensure right number of bins
	if (_data.size() != _waveset.nTbin()){
		_data = std::vector<std::vector<std::vector<double> > >(_waveset.nTbin());
	};
	if (_data[tbin].size() != _waveset.nBins()){
		_data[tbin] = std::vector<std::vector<double> >(_waveset.nBins());
	};
	_data[tbin][bin] = data;
	update_is_active();

	if (data.size() == _waveset.nPoints()*_waveset.nPoints()){
		return true;
	};
	return false;
};
//########################################################################################################################################################
///Set coma manually
bool full_covariance::set_coma(
							int 								tbin,
							int 								bin,
							std::vector<std::vector<double> > 				coma){

	//print_matrix(coma);
	// Ensure right number of bins
	if(_coma.size() != _waveset.nTbin()){
		_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_waveset.nTbin());
	};
	if(_coma[tbin].size() != _waveset.nBins()){
		_coma[tbin] = std::vector<std::vector<std::vector<double> > >(_waveset.nBins());
	};
	_coma[tbin][bin] = coma;
	if (coma.size() == _waveset.nPoints()*_waveset.nPoints()){
		for (size_t i=0;i<coma.size();i++){
			if (coma[i].size() != coma.size()){
				std::cout<<"Warning: Set coma is not quadratic."<<std::endl;
				return false;
			};
		};
		return true;
	};
	return false;
};
//########################################################################################################################################################
///Loads the data for a t' bin from a file
void full_covariance::loadData(
							int 								tbin,
							const char* 							dataFile){

	_data[tbin] = std::vector<std::vector<double> >();
	std::fstream data(dataFile,std::ios_base::in);
	double val;
	size_t nDat = _waveset.nPoints()*_waveset.nPoints();
	std::vector<double> data_bin;
	while(data >> val){
		data_bin.push_back(val);
		if (data_bin.size() == nDat){
			_data[tbin].push_back(data_bin);
			data_bin = std::vector<double>();
		};
	};
	update_is_active();
	if (_waveset.nBins() != _data[tbin].size()){
		std::cout << "Warning: _waveset.nBins()="<<_waveset.nBins()<<" != _data.size()="<<_data[tbin].size()<<std::endl;
	}else{
		std::cout << "File delivered the right size for _data"<<std::endl;
	};
};
//########################################################################################################################################################
///Loads the (inverted) covariance matrix for a t' bin from 'comaFile'
void full_covariance::loadComa(
							int 								tbin,
							const char* 							comaFile){

	_coma[tbin]=std::vector<std::vector<std::vector<double> > >();
	std::fstream data(comaFile,std::ios_base::in);
	double val;
	size_t nDat = _waveset.nPoints()*_waveset.nPoints();
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
	if (_waveset.nBins() != _coma[tbin].size()){
		std::cout << "Warning: _waveset.nBins()="<<_waveset.nBins()<<" != _coma.size()="<<_coma[tbin].size() << std::endl;
	}else{
		std::cout<< "File delivered the right size for _coma"<<std::endl;
	};
};
//########################################################################################################################################################
///Returns the number of couplings
int full_covariance::getNcpl(){
	std::vector<size_t> nnn = *_waveset.n_cpls();
	size_t n = nnn.size();
	size_t maxx = 0;
	for (size_t i=0;i<n;++i){
		if (nnn[i] > maxx){
			maxx = nnn[i];
		};
	};
	maxx+=1;
	maxx *= _waveset.nTbin();
	return maxx;
};
//########################################################################################################################################################
/// Sets br and cpl in one t' bin so, that it corresponds tothe full cpls
std::vector<std::vector<std::complex<double> > > full_covariance::full_to_br_cpl(std::vector<std::complex<double> > &cpl){

	size_t nFtw  = _waveset.nFtw();
	size_t nTbin = _waveset.nTbin();
	size_t nCpl  = _nCpl;
	size_t cpls  = nCpl/nTbin;
	size_t nBra  = _waveset.nBranch();

	std::vector<std::vector<std::complex<double> > > ret;
	ret.push_back(std::vector<std::complex<double> >(cpls,std::complex<double>(0.,0.)));
	ret.push_back(std::vector<std::complex<double> >(nBra,std::complex<double>(0.,0.)));

	std::vector<int> n_branch= *_waveset.n_branch();
	std::vector<size_t> n_cpls = *_waveset.n_cpls();

	std::complex<double> phase = conj(cpl[0])/abs(cpl[0]);

	for (size_t func=0;func<nFtw;++func){
		int nBr = n_branch[func];
		size_t nCp = n_cpls[func];
		if (nBr ==-1){
			ret[0][nCp] = cpl[func]*phase;
		}else{
			if(ret[0][nCp] == std::complex<double>(0.,0.)){
				ret[0][nCp] = cpl[func]*phase;
				ret[1][nBr] = std::complex<double>(1.,0.);
			}else{
				ret[1][nBr] = cpl[func]*phase/ret[0][nCp];
			};
		};
	};
	return ret;
};
//########################################################################################################################################################
///Gets the total number of paramters with only anchor couplings
int full_covariance::getNtot(){

	return 2*getNcpl() + _waveset.getNpar() + 2*_waveset.nBranch() + _waveset.getNiso();
};
//########################################################################################################################################################
///Sets all data to 0.
void full_covariance::nullify(){

	for (size_t tbin = 0; tbin < _waveset.nTbin();tbin++){
		for (size_t bin =0; bin<_waveset.nBins();bin++){
			int nDat = _waveset.nPoints()*_waveset.nPoints();
			for (int i=0;i<nDat;i++){
				_data[tbin][bin][i]=0.;
			};
		};
	};
};
//########################################################################################################################################################
///Updates the number of couplings
void full_covariance::update_n_cpls(){

	std::vector<size_t> nnn = *_waveset.n_cpls();
	size_t n = nnn.size();
	size_t maxx = 0;
	for (size_t i=0;i<n;++i){
		if (nnn[i] > maxx){
			maxx = nnn[i];
		};
	};
	maxx+=1;
	_nBrCpl =  maxx;
};
//########################################################################################################################################################
///Updates the internal ranges for evaluation
void full_covariance::update_min_max_bin(){

// Override the mthod in waveset, because the limits are given by the anchor wave in this method
	_waveset.setMinBin(0);
//	std::cout<<"Call update_min_max_bin()"<<std::endl;
	_waveset.setMaxBin((*_waveset.binning()).size());
	if ((*_waveset.lowerLims()).size() == 0 or (*_waveset.upperLims()).size() == 0){
		std::cout<< "Warning: Wave limits not set, omitting internal setting of _waveset.minBin() and _waveset.maxBin()"<<std::endl;
		std::cout<< "Warning: Setting _waveset.minBin() = 0; _waveset.maxBin() = _binning.size() = "<<(*_waveset.binning()).size()<<std::endl;
		_waveset.setMinBin(0);
		_waveset.setMaxBin((*_waveset.binning()).size());
		return;
	};
	double ancMin = (*_waveset.lowerLims())[0];
	double ancMax = (*_waveset.upperLims())[0];
//	std::cout<<ancMin<<"-"<<ancMax<<std::endl;
	if ((*_waveset.binning()).size() == 0){
		std::cout<< "Warning: No binning set, cannot determine _waveset.minBin(), _waveset.maxBin()" <<std::endl;
	};
	for (size_t i=0;i<(*_waveset.binning()).size()-1;i++){
		double up = (*_waveset.binning())[i+1];
		double low= (*_waveset.binning())[i];
		if (ancMin >= low and ancMin <up){
			_waveset.setMinBin(i);
		};
		if (ancMax > low and ancMax <=up){
			_waveset.setMaxBin(i+1);
		};
	};
	update_n_cpls();
//	std::cout<<"_waveset.minBin(): "<<_waveset.minBin()<<std::endl;
//	std::cout<<"_waveset.maxBin(): "<<_waveset.maxBin()<<std::endl;
};
//########################################################################################################################################################
/// Prints the internal status
void full_covariance::printStatus()																const{
	_waveset.printStatus();
	std::cout<<std::endl<<"FULL COVARIANCE:"<<std::endl;
	std::cout<<std::endl<<"PARAMETER NUMBERS:"<<std::endl;
	std::cout<<"_nTot: "<<_nTot<<"; _nPar: "<<_nPar<<"; _nCpl: "<<_nCpl<<std::endl;
	std::cout<<"_nBra: "<<_nBra<<"; _nIso: "<<_nIso<< "; _nBrCpl: "<<_nBrCpl<<std::endl;
	std::cout<<std::endl<<std::endl<<"PARAMETERS:"<<std::endl<<"_parameters"<<std::endl;
	print_vector(_parameters);
	std::cout<<std::endl<<"_upper_parameter_limits"<<std::endl;
	print_vector(_upper_parameter_limits);
	std::cout<<std::endl<<"_lower_parameter_limits"<<std::endl;
	print_vector(_lower_parameter_limits);
	std::cout<<std::endl<<"_parNames"<<std::endl;
	print_vector(_parNames);
	std::cout<<std::endl<<"Omit printing of _data and _coma"<<std::endl;
	std::cout<<std::endl<<"OTHER MEMBERS:"<<std::endl;
	std::cout<<std::endl<<"_nOut: "<<_nOut<<std::endl;
	std::cout<<std::endl<<"_count: "<<_count<<std::endl;

};
//########################################################################################################################################################
///Complex conjugates _data, changes _coma accordingly
void full_covariance::conjugate(){ //<<need>>

	size_t nPoints = _waveset.nPoints();
	for (size_t ii = 0;ii<nPoints*nPoints;++ii){
		size_t i = ii/nPoints;
		size_t j = ii-nPoints*i;
		if (j<i){
			for(size_t tbin=0;tbin<_waveset.nTbin();++tbin){
				for(size_t bin=0;bin<_waveset.nBins();++bin){
//					_data[tbin][bin][ii]*=-1;
					for (size_t jj=0;jj<nPoints*nPoints;++jj){
						_coma[tbin][bin][ii][jj]*=-1;
						_coma[tbin][bin][jj][ii]*=-1;
					};
				};
			};
		};
	};
	for (size_t i=0;i<nPoints;++i){
		for (size_t j=0;j<i;++j){
			for(size_t tbin=0;tbin<_waveset.nTbin();++tbin){
				for(size_t bin=0;bin<_waveset.nBins();++bin){
					_data[tbin][bin][i*nPoints+j]*=-1;
				};
			};
		};
	};

};
//########################################################################################################################################################
#ifdef USE_YAML
///Loads _data and _coma from a YAML file
bool full_covariance::loadDataComa(
							YAML::Node					&waveset){


	if (waveset["data_files"].size() != _waveset.nTbin() ){
		std::cerr<<"Error: Number of data-files does not match number of t' bins"<<std::endl;
	};
	if (waveset["coma_files"].size() != _waveset.nTbin()){
		std::cerr<<"Error: Number of coma-files does not match number of t' bins"<<std::endl;
	};
	update_min_max_bin(); //Has to be called to set internal definitions right, since not all belong to the same class now
	update_definitions(); // Has to be called once with all functions set, to get internal handling of parameter numbers right
	if (waveset["data_files"].size() == _waveset.nTbin() and waveset["coma_files"].size() == _waveset.nTbin()){
		for(size_t i=0;i<_waveset.nTbin();i++){
			loadData(i,waveset["data_files"][i].as<std::string>().c_str());
			loadComa(i,waveset["coma_files"][i].as<std::string>().c_str());
		};
	}else{
		std::cerr<<"Error: Wrong number of files, no _data or _coma loaded"<<std::endl;
		return false;
	};
	if (waveset["conjugate_fit_result"]){
		if (waveset["conjugate_fit_result"].as<bool>()){
			conjugate();
		};
	};
	return true;
};
//########################################################################################################################################################
///Loads paramters values and limits from a YAML file
bool full_covariance::loadParameterValues(
							YAML::Node					&waveset,
							YAML::Node					&param){

	std::map<std::string,int> fMap;
	int fCount = 1; // Starts as one, since map[nonExistingKey] = 0 -> use this to determine, whether a function exists.
	int iCount = 1;
	int pCount = 0;
	int ipCount = 0;
	int nWaves = waveset["waves"].size();
	bool ookk = true;
	for (int i=0;i<nWaves;i++){
		int nFunc = waveset["waves"][i]["parametrizations"].size();
		for (int j=0;j<nFunc;j++){
			std::string fName;
			std::string iName;
			bool isisobar = false;
			if (waveset["waves"][i]["parametrizations"][j].size()==0){
				fName = waveset["waves"][i]["parametrizations"][j].as<std::string>();
//				std::cout<<fName<<" alone"<<std::endl;
			}else if(waveset["waves"][i]["parametrizations"][j].size()==2) {
				isisobar = true;
				fName = waveset["waves"][i]["parametrizations"][j][0].as<std::string>();
				iName = waveset["waves"][i]["parametrizations"][j][1].as<std::string>();
//				std::cout<<iName<<" to "<<fName<<std::endl;
			};
			if(param[fName]){
				if (fMap[fName]==0){
					fMap[fName]=fCount;
					fCount++;
					int nPar = param[fName]["parameters"].size();
					for(int par=0;par<nPar;par++){
						double value = param[fName]["parameters"][par]["value"].as<double>();
						if(param[fName]["parameters"][par]["upper_limit"] and param[fName]["parameters"][par]["lower_limit"]){
							double upper = param[fName]["parameters"][par]["upper_limit"].as<double>();
							double lower = param[fName]["parameters"][par]["lower_limit"].as<double>();
							setParLimits(2*_nCpl+pCount,upper,lower);
						};
						setParameter(2*_nCpl+pCount,value);
						pCount++;
					};
				};
			}else{
				ookk = false;
				std::cerr << "Error: '"<<fName<<"' not defined in parametrization file, no values set"<<std::endl;
			};
			if (isisobar){
				if(param[iName]){
					if(fMap[iName]==0){
						fMap[iName]=iCount;
						iCount++;
						int nPar =  param[iName]["parameters"].size();
						for(int par=0;par<nPar;par++){
							double value = param[iName]["parameters"][par]["value"].as<double>();
							if(param[iName]["parameters"][par]["upper_limit"] and param[iName]["parameters"][par]["lower_limit"]){
								double upper = param[iName]["parameters"][par]["upper_limit"].as<double>();
								double lower = param[iName]["parameters"][par]["lower_limit"].as<double>();
								setParLimits(2*_nCpl+_nPar+2*_nBra+ipCount,upper,lower);
							};
							setParameter(2*_nCpl+_nPar+2*_nBra+ipCount,value);
							ipCount++;
						};
					};
				}else{
					ookk = false;
				};
			};
		};
	};
	// Load explicit definitions
	if (waveset["parameters"]){
		int nDef = waveset["parameters"].size();
		for (int i=0;i<nDef;i++){
			std::string name = waveset["parameters"][i]["name"].as<std::string>();
			int par = getParNumber(name);
			if (par==-1){
				std::cout<<"Warning: Unknown parameter in YAML file: "<<name<<std::endl;
				ookk = false;
			}else{
				if (waveset["parameters"][i]["value"]){
					setParameter(name,waveset["parameters"][i]["value"].as<double>());
				};
				if (waveset["parameters"][i]["upper_limit"] and waveset["parameters"][i]["lower_limit"]){
					setParLimits(par,waveset["parameters"][i]["upper_limit"].as<double>(),waveset["parameters"][i]["lower_limit"].as<double>());
				};
			};
		};
	};
	return ookk;
};
#endif//USE_YAML
//########################################################################################################################################################
///Initializes _data and _coma, when setting the t' binning
void full_covariance::setTbinning(std::vector<std::vector<double> > binning){
	_waveset.setTbinning(binning);
	if (_data.size() != _waveset.nTbin()){
		if (_data.size() != 0){
			std::cout<<"Warning: _data.size() != _waveset.nTbin(), but nonzero. Previous set _data[][][] will be lost."<<std::endl;
		};
		_data = std::vector<std::vector<std::vector<double> > >(_waveset.nTbin());
	};
	if (_coma.size() != _waveset.nTbin()){
		if (_coma.size() != 0){
			std::cout<<"Warning: _coma.size() != _waveset.nTbin(), bun nonzero. Previous set _coma[][][][] will be lost."<<std::endl;
		};
		_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_waveset.nTbin());
	};
};
//########################################################################################################################################################
///Updates internal definitions
void full_covariance::update_definitions(){

	_nTot = getNtot();
	_nPar = _waveset.getNpar();
	_nCpl = getNcpl();
	_nBra = _waveset.nBranch();
	_nIso = _waveset.getNiso();
	std::vector<std::string> names;
	std::stringstream str_count;
	for (size_t i=0;i<_nCpl;i++){
		str_count<<i;
		names.push_back("reC"+str_count.str());
		names.push_back("imC"+str_count.str());
		str_count.str("");
	};
	for (size_t i=0;i<_nPar;i++){
		names.push_back(_waveset.getParameterName(i));
	};
	for (size_t i=0;i<_nBra;i++){
		str_count<<i;
		names.push_back("reB"+str_count.str());
		names.push_back("imB"+str_count.str());
		str_count.str("");
	};
	for (size_t i=0;i<_nIso;i++){
		names.push_back(_waveset.getIsoParName(i));
	};
	_parNames = names;
};
//########################################################################################################################################################
///Updates, which data-point is actually active
void full_covariance::update_is_active(){ //<<need>>

	_is_active = std::vector<std::vector<bool> >(_waveset.nBins(),std::vector<bool>(_waveset.nPoints(),true));
	size_t  nPoint = _waveset.nPoints();
	std::vector<int> point_to_wave = *(_waveset.point_to_wave());
	for(size_t bin=0;bin<_waveset.nBins();++bin){
		double mass = _waveset.get_m(bin);
		for (size_t point = 0;point<nPoint;++point){
			size_t wave = point_to_wave[point];
			double upperLim = (*_waveset.upperLims())[wave];
			double lowerLim = (*_waveset.lowerLims())[wave];
			if (mass > upperLim or mass < lowerLim){
				_is_active[bin][point] = false;
			};
		};
	};
};
//########################################################################################################################################################
///Writes plots with the internal _parameters
void full_covariance::write_plots(
							std::string						filename,
							int							tbin)					const{

	write_plots(filename,tbin,parameters());
};
//########################################################################################################################################################
///Write plots with the usual parameter ordering
void full_covariance::write_plots(
							std::string						filename,
							int							tbin,
							const std::vector<double>				&param)					const{

	std::vector<std::complex<double> > cpl(_nCpl);
	std::vector<double> par(_nPar);	
	std::vector<std::complex<double> > bra(_nBra);
	std::vector<double> iso(_nIso);
	int count =0;
	for (size_t i=0;i<_nCpl;i++){
		cpl[i] = std::complex<double>(param[count],param[count+1]);
		count+=2;
	};
	for (size_t i=0;i<_nPar;i++){
		par[i] = param[count];
		count++;
	};
	for (size_t i=0;i<_nBra;i++){
		bra[i] = std::complex<double>(param[count],param[count+1]);
		count+=2;
	};
	for (size_t i=0;i<_nIso;i++){
		iso[i]=param[count];
		count++;
	};
	write_plots(filename,tbin,cpl,par,bra,iso);
};
//########################################################################################################################################################
///Writes plots for one t' bin, only diagonal elements at the moment
void full_covariance::write_plots(
							std::string						filename,
							int 							tbin,
							const std::vector<std::complex<double> >		&cpl,
							const std::vector<double>				&par,
							const std::vector<std::complex<double> >		&bra,
							const std::vector<double> 				&iso)					const{

	std::cout<<std::endl;
	std::vector<std::complex<double> > cpl_all = getAllCouplings(tbin,cpl,par,bra,iso);
	std::cout<<std::endl;
	std::vector<std::vector<std::complex<double> > > iso_eval;
	if(_waveset.has_isobars()){
		iso_eval = _waveset.iso_funcs(&iso[0]);
	};
	std::ofstream write_out;
	write_out.open(filename.c_str());
	std::cout<<"write_plots(...): Chi2 for the used paramters is: "<<EvalTbin(tbin,&cpl_all[0],&par[0],&iso[0])<<std::endl;//[0]//
	for (size_t bin=0;bin<_waveset.nBins();bin++){
		double mass = _waveset.get_m(bin);
		std::vector<double> var = _waveset.getVar(mass,tbin);
		std::vector<std::complex<double> > amplitudes = _waveset.amps(&var[0],&cpl_all[0],&par[0],iso_eval);
		write_out<<mass<<" ";
		for (size_t i=0; i<_waveset.nPoints(); i++){
			for (size_t j=0; j<_waveset.nPoints(); j++){
				std::complex<double> inter = amplitudes[i]*conj(amplitudes[j]);
				if (i==j){
					write_out<<inter.real()<<" "<<_data[tbin][bin][(_waveset.nPoints()+1)*i]<<" "<<pow(_coma[tbin][bin][(_waveset.nPoints()+1)*i][(_waveset.nPoints()+1)*i],-.5)<<" ";
				};
				if (i<j){
					write_out<<inter.real()<<" "<<_data[tbin][bin][_waveset.nPoints()*i+j]<<" "<<pow(_coma[tbin][bin][_waveset.nPoints()*i+j][_waveset.nPoints()*i+j],-.5)<<" ";
				};
				if(i>j){
					write_out<<inter.imag()<<" "<<_data[tbin][bin][_waveset.nPoints()*i+j]<<" "<<pow(_coma[tbin][bin][_waveset.nPoints()*i+j][_waveset.nPoints()*i+j],-.5)<<" ";
				};
			};
		};
		write_out<<std::endl;
	};
	write_out.close();
};
//########################################################################################################################################################
