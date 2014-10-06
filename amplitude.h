#ifndef BREIT_WIGNERS_WEIT_BRIGNERS
#define BREIT_WIGNERS_WEIT_BRIGNERS
#include<string>
#include<complex>
#include<iostream>
#include"amplitude_functions.h"

class amplitude {
	public:
		amplitude();
		amplitude(size_t nVar, size_t nPar, size_t nCon, int funcId);

		virtual std::string 	type()									const		{return "constant_function";};
		std::string		name()									const		{return _name;};

		std::complex<double> Eval(const double* var)							const		{return Eval(var, _parameters, _constants);};
		std::complex<double> Eval(const double* var, const double* par)					const		{return Eval(var, par, _constants);};
		virtual std::complex<double> Eval(const double* var, const double* par, const double* con)	const		{return std::complex<double>(1.,0.);};

#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var)							const		{return Eval(var, (adouble*)_parameters, (adouble*)_constants);};
		std::complex<adouble> Eval(const adouble* var, const adouble* par)				const		{return Eval(var, par, (adouble*)_constants);};
		virtual std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)	const		{return std::complex<adouble>(1.,0.);};
#endif//ADOL_ON

		size_t 			nVar()									const		{return _nVar;};
		size_t 			nPar()									const		{return _nPar;};
		size_t 			nCon()									const		{return _nCon;};

		int			L()									const		{return _L;};
		void			setL(int L)										{_L=L;};

		bool 			setPar(size_t n, double val);
		bool 			setCon(size_t n, double val);

		bool			setPars(const double* vals);
		bool			setCons(const double* vals);

		const double*		parameters()								const		{return _parameters;};
		const double*		constants()								const		{return _constants;};
	
		double			getParameter(size_t n)							const;
		double			getConstant(size_t n)							const;

		std::string		getVarType(size_t n)							const;					
		std::string		getParType(size_t n)							const;
		std::string		getConType(size_t n)							const;
		std::string		getParName(size_t n)							const;
		std::string		getConName(size_t n)							const;

		int			getFuncId()								const		{return _funcId;};

		void			setFunctionName(std::string name)							{_name = name;};
		bool			setParName(size_t n, std::string name);
		bool			setConName(size_t n, std::string name);

	protected:

		int 			_L;		// Spin, treat special here, since its essential

		std::string*		_var_types;	// Type of the variables // Convention: {m, t',...}
		std::string*		_par_types;	// Names of the parameters
		std::string*		_con_types;	// Names of the constants

		std::string		_name;		// Name of the functions
		std::string*		_par_names;	// Names of the parameters
		std::string*		_con_names;	// Names of the constants

		const size_t		_nVar;		// # of Variables
		const size_t		_nPar;		// # of Paramters
		const size_t		_nCon;		// # of Constants

		const int		_funcId;	// Id # of the function type, to be cosistent with the old method

		double*			_parameters;	// Paramters
		double*			_constants;	// Constants

};

amplitude::amplitude():_nVar(0),_nPar(0),_nCon(0),_funcId(-1){ // Id for constant function is -1

		_name = "unnamed_constant_function";
};

amplitude::amplitude(size_t nVar, size_t nPar, size_t nCon, int funcId):_nVar(nVar), _nPar(nPar), _nCon(nCon), _funcId(funcId){

		_name = "unnamed_constant_function";
};

bool amplitude::setPar(size_t n, double val){
	if (n < _nPar){
		_parameters[n] = val;
		return true;
	}else{
		std::cerr << "Error: Can't set parameter #"<<n<<" for "<<type()<<std::endl;		
		return false;
	};
};

bool amplitude::setCon(size_t n, double val){
	if (n < _nCon){
		_constants[n] = val;
		return true;
	}else{
		std::cerr << "Error: Can't set constant #"<<n<<" for "<<type()<<std::endl;		
		return false;
	};
};

bool amplitude::setPars(const double* vals){
	bool ret = true;
	for (size_t iii =0;iii<_nPar;iii++){
		if (not setPar(iii,vals[iii])){
			ret = false;
		};
	};
	return ret;
};

bool amplitude::setCons(const double* vals){
	bool ret = true;
	for (size_t iii =0;iii<_nCon;iii++){
		if (not setCon(iii,vals[iii])){
			ret = false;
		};
	};
	return ret;
};

std::string amplitude::getVarType(size_t n)		const{
	if (n<_nVar){
		return _var_types[n];
	}else{
		std::cerr <<"Error: getVarType(...): _nVar  = "<<_nVar<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
std::string amplitude::getParType(size_t n)		const{
	if (n<_nPar){
		return _par_types[n];
	}else{
		std::cerr <<"Error: getParType(...): _nPar  = "<<_nPar<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
std::string amplitude::getConType(size_t n)		const{
	if (n<_nCon){
		return _con_types[n];
	}else{
		std::cerr <<"Error: getConType(...): _nCon  = "<<_nCon<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
std::string amplitude::getParName(size_t n)		const{
	if (n<_nPar){
		return _par_names[n];
	}else{
		std::cerr <<"Error: getParName(...): _nPar  = "<<_nPar<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
std::string amplitude::getConName(size_t n)		const{
	if (n<_nCon){
		return _con_names[n];
	}else{
		std::cerr <<"Error: getConName(...): _nCon  = "<<_nCon<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
bool amplitude::setParName(size_t n, std::string name){
	if(n<_nPar){
		_par_names[n] = name;
		return true;
	}else{
		std::cerr<<"Error: Can't set _par_names["<<n<<"] for "<<type()<<std::endl;
		return false;
	};
};
bool amplitude::setConName(size_t n, std::string name){
	if(n<_nCon){
		_con_names[n] = name;
		return true;
	}else{
		std::cerr<<"Error: Can't set _con_names["<<n<<"] for "<<type()<<std::endl;
		return false;
	};
};

////////	////////	END OF BASE CLASS DEFINITION	////////	BEGIN SPECIAL DEFINITIONS	/////////// 

////////	////////	constant_function 		///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef amplitude constant_function; // Amplitude base class is already the constant function


////////	////////	breit_wigner 			///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class breit_wigner : public amplitude{

	public:
		breit_wigner();

		std::string type()											const		{return "breit_wigner";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

breit_wigner::breit_wigner():amplitude(1,2,0,0){

	_name = "unnamed_breit_wigner";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[2] = {"mass","width"};
	_par_types = par_types;

	std::string par_names[2] = {"unnamed_mass","unnamed_width"};
	_par_names = par_names;
};

template <typename xdouble>
std::complex<xdouble> breit_wigner::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const{
	std::complex<xdouble> denominator = std::complex<xdouble>(par[0]*par[0]-var[0]*var[0],-par[0]*par[1]);
	return std::complex<xdouble>(par[0]*par[1])/denominator;
};

////////	////////	mass_dep_breit_wigner		///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class mass_dep_breit_wigner : public amplitude{

	public:
		mass_dep_breit_wigner();

		std::string type()											const		{return "mass_dep_breit_wigner";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

mass_dep_breit_wigner::mass_dep_breit_wigner():amplitude(1,2,2,1){

	_name = "unnamed_mass_dep_breit_wigner";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[2] = {"mass","width"};
	_par_types = par_types;

	std::string par_names[2] = {"unnamed_mass","unnamed_width"};
	_par_names = par_names;

	std::string con_types[2] = {"m_Pi","m_Iso"};
	_con_types = con_types;

	std::string con_names[2] = {"m_Pi","m_Iso"};
	_con_names=con_names;
};

template <typename xdouble>
std::complex<xdouble> mass_dep_breit_wigner::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)	const{
		xdouble m   = var[0];

		xdouble m0  = par[0];
		xdouble G0  = par[1];

		xdouble mPi = con[0];
		xdouble mIso= con[1];

		xdouble q0 = breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso*mIso);
		xdouble q  = breakupMomentumReal<xdouble>(m* m ,mPi*mPi,mIso*mIso);
		xdouble Fl = barrierFactor<xdouble>(q,_L);
		xdouble Fl0= barrierFactor<xdouble>(q0,_L);

		xdouble G  = G0* m0/m * q*Fl*Fl/q0/Fl0/Fl0; //G0 * m0/m q*Fl^2/(q0*Fl0^2)
		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);
		return std::complex<xdouble>(m0*G0,0.)/denominator;	
};

////////	////////	two_channel_breit_wigner	///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class two_channel_breit_wigner : public amplitude{

	public:
		two_channel_breit_wigner();

		std::string type()											const		{return "two_channel_breit_wigner";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

two_channel_breit_wigner::two_channel_breit_wigner():amplitude(1,2,4,2){

	_name = "unnamed_two_channel_breit_wigner";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[2] = {"mass","width"};
	_par_types = par_types;

	std::string par_names[2] = {"unnamed_mass","unnamed_width"};
	_par_names = par_names;

	std::string con_types[4] = {"m_Pi","m_Iso1","m_Iso2","branching"};
	_con_types = con_types;

	std::string con_names[4] = {"m_Pi","m_Iso1","m_Iso2","branching"};
	_con_names=con_names;
};

template <typename xdouble>
std::complex<xdouble> two_channel_breit_wigner::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];

		xdouble m0    = par[0];
		xdouble G0    = par[1];

		xdouble mPi   = con[0];
		xdouble mIso1 = con[1];
		xdouble mIso2 = con[2];
		xdouble X  = con[3];

		xdouble R = 5.;

		xdouble psl1 = psl<xdouble>(m, mPi, mIso1, R, _L);
		xdouble psl2 = psl<xdouble>(m, mPi, mIso2, R, _L);

		xdouble psl10= psl<xdouble>(m0, mPi, mIso1, R, _L);
		xdouble psl20= psl<xdouble>(m0, mPi, mIso2, R, _L);

		xdouble G = G0 * m0/m * ((1-X) * psl1/psl10 + X * psl2/psl20);

		std::complex<xdouble> denominator  = std::complex<xdouble>(m0*m0-m*m,-m0*G);
		return std::complex<xdouble>(m0*G0,0)/denominator;
};

////////	////////	vandermeulen_phase_space	///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class vandermeulen_phase_space : public amplitude{

	public:
		vandermeulen_phase_space();

		std::string type()											const		{return "vandermeulen_phase_space";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

vandermeulen_phase_space::vandermeulen_phase_space():amplitude(1,1,2,3){

	_name = "unnamed_vandermeulen_phase_space";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[1] = {"alpha"};
	_par_types = par_types;

	std::string par_names[1] = {"unnamed_alpha"};
	_par_names = par_names;

	std::string con_types[2] = {"m_Pi","m_Iso"};
	_con_types = con_types;

	std::string con_names[2] = {"m_Pi","m_Iso"};
	_con_names=con_names;
};

template <typename xdouble>
std::complex<xdouble> vandermeulen_phase_space::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];

		xdouble alpha = par[0];

		xdouble mPi   = con[0];
		xdouble mIso  = con[1];

		xdouble ampor = mPi + mIso;
		std::complex<xdouble> value;

		if ( m > ampor){
			xdouble S = m*m;
			xdouble E = (S + mPi * mPi - mIso*mIso)/(2*m);
			xdouble PSQ = E*E - mPi*mPi;
			value = std::complex<xdouble>(exp(alpha*PSQ),0.);
		}else{
			value = std::complex<xdouble>(1.,0.);			
		};
		return value;
};

////////	////////	valera_dorofeev_background	///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class valera_dorofeev_background : public amplitude{

	public:
		valera_dorofeev_background();

		std::string type()											const		{return "valera_dorofeev_background";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

valera_dorofeev_background::valera_dorofeev_background():amplitude(1,2,1,4){

	_name = "unnamed_valera_dorofeev_background";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[2] = {"alpha","beta"};
	_par_types = par_types;

	std::string par_names[2] = {"unnamed_alpha","unnamed_beta"};
	_par_names = par_names;

	std::string con_types[1] = {"m_0"};
	_con_types = con_types;

	std::string con_names[1] = {"m_0"};
	_con_names=con_names;
};

template <typename xdouble>
std::complex<xdouble> valera_dorofeev_background::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];

		xdouble alpha = par[0];
		xdouble beta  = par[1];

		xdouble m0    = con[0];		
		return std::complex<xdouble>(pow((m-m0)/0.5,alpha)*exp(-beta*(m-m0-0.5)),0);
};

////////	////////	bowler_parametrization		///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class bowler_parametrization : public amplitude{

	public:
		bowler_parametrization();

		std::string type()											const		{return "bowler_parametrization";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

bowler_parametrization::bowler_parametrization():amplitude(1,2,0,5){

	_name = "unnamed_bowler_parametrization";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[2] = {"mass","width"};
	_par_types = par_types;

	std::string par_names[2] = {"unnamed_mass","unnamed_width"};
	_par_names = par_names;
};

template <typename xdouble>
std::complex<xdouble> bowler_parametrization::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];

		xdouble m0    = par[0];
		xdouble G0    = par[1];

		xdouble G = G0*	bowler_integral_table<xdouble>(m)/bowler_integral_table<xdouble>(m0) * m0/m;
		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);

		return std::complex<xdouble>(sqrt(m0*G0),0)/denominator;
};
////////	////////	flatte				///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class flatte : public amplitude{

	public:
		flatte();

		std::string type()											const		{return "flatte";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

flatte::flatte():amplitude(1,3,2,6){

	_name = "unnamed_flatte";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[3] = {"mass","g_1","g_2"};
	_par_types = par_types;

	std::string par_names[3] = {"unnamed_mass","unnamed_g_1","unnamed_g_2"};
	_par_names = par_names;

	std::string con_types[2] = {"m_Pi","m_K"};
	_con_types = con_types;

	std::string con_names[2] = {"m_Pi","m_K"};
	_con_names=con_names;
};

template <typename xdouble>
std::complex<xdouble> flatte::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];

		xdouble m0    = par[0];
		xdouble g1    = par[1];
		xdouble g2    = par[2];

		xdouble mPi   = con[0];
		xdouble mK    = con[1];
		
		xdouble qpp= breakupMomentumReal<xdouble>(m*m,mPi*mPi,mPi*mPi);
		xdouble qKK= breakupMomentumReal<xdouble>(m*m,mK*mK,mK*mK);

		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-(g1*qpp*qpp + g2*qKK*qKK));
		return std::complex<xdouble>(1,0)/denominator;
};

////////	////////	gaus				///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class gaus : public amplitude{

	public:
		gaus();

		std::string type()											const		{return "gaus";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

gaus::gaus():amplitude(1,2,0,9){

	_name = "unnamed_gaus";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[2] = {"mu","sigma"};
	_par_types = par_types;

	std::string par_names[2] = {"unnamed_mu","unnamed_sigma"};
	_par_names = par_names;

};

template <typename xdouble>
std::complex<xdouble> gaus::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];

		xdouble m0    = par[0];
		xdouble sig   = par[1];

		return std::complex<xdouble>(exp(-(m-m0)*(m-m0)/2/sig/sig),0);
};
////////	////////	polynomial			///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class polynomial : public amplitude{

	public:
		polynomial();

		std::string type()											const		{return "polynomial";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

polynomial::polynomial():amplitude(1,5,0,10){

	_name = "unnamed_polynomial";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[5] = {"c0","c1","c2","c3","c4"};
	_par_types = par_types;

	std::string par_names[5] = {"unnamed_c0","unnamed_c1","unnamed_c2","unnamed_c3","unnamed_c4"};
	_par_names = par_names;

};

template <typename xdouble>
std::complex<xdouble> polynomial::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];

		xdouble c0    = par[0];
		xdouble c1    = par[1];
		xdouble c2    = par[2];
		xdouble c3    = par[3];
		xdouble c4    = par[4];

		xdouble ret = c4*m*m*m*m + c3*m*m*m + c2*m*m + c1*m + c0;
		return std::complex<xdouble>(ret,0);
};

////////	////////	mass_dep_bw_2			///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class mass_dep_bw_2 : public amplitude{

	public:
		mass_dep_bw_2();

		std::string type()											const		{return "mass_dep_bw_2";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

mass_dep_bw_2::mass_dep_bw_2():amplitude(1,2,4,22){

	_name = "unnamed_mass_dep_bw_2";

	std::string var_types[1] = {"mass"};
	_var_types = var_types;

	std::string par_types[2] = {"mass","width"};
	_par_types = par_types;

	std::string par_names[2] = {"unnamed_mass","unnamed_width"};
	_par_names = par_names;

	std::string con_types[4] = {"m_Pi","m_Iso1","m_Iso2","branching"};
	_con_types = con_types;

	std::string con_names[4] = {"m_Pi","m_Iso1","m_Iso2","branching"};
	_con_names=con_names;
};

template <typename xdouble>
std::complex<xdouble> mass_dep_bw_2::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];

		xdouble m0    = par[0];
		xdouble G0    = par[1];

		xdouble mPi   = con[0];
		xdouble mIso1 = con[1];
		xdouble mIso2 = con[2];
		xdouble X     = con[3];

		xdouble q1 = breakupMomentumReal<xdouble>(m*m,mPi*mPi,mIso1*mIso1);
		xdouble q10= breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso1*mIso1);
		xdouble q2 = breakupMomentumReal<xdouble>(m*m,mPi*mPi,mIso2*mIso2);
		xdouble q20= breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso2*mIso2);
		xdouble Fl1= barrierFactor<xdouble>(q1,_L);
		xdouble Fl10=barrierFactor<xdouble>(q10,_L);
		xdouble Fl2= barrierFactor<xdouble>(q2,_L);
		xdouble Fl20=barrierFactor<xdouble>(q20,_L);	

		xdouble G  = G0 * m0/m* ((1-X)*q1*Fl1*Fl1/q10/Fl10/Fl10 + X* q2*Fl2*Fl2/q20/Fl20/Fl20);
		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);
		return std::complex<xdouble>(m0*G0,0.)/denominator;	
};

////////	////////	t_dependent_background		///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class t_dependent_background : public amplitude{

	public:
		t_dependent_background();

		std::string type()											const		{return "t_dependent_background";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const xdouble* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const adouble* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

t_dependent_background::t_dependent_background():amplitude(2,4,3,101){

	_name = "unnamed_t_dependent_background";

	std::string var_types[2] = {"mass","t'"};
	_var_types = var_types;

	std::string par_types[4] = {"b","c0","c1","c2"};
	_par_types = par_types;

	std::string par_names[4] = {"unnamed_b","unnamed_c0","unnamed_c1","unnamed_c2"};
	_par_names = par_names;

	std::string con_types[3] = {"m_0","m_Pi","mIso"};
	_con_types = con_types;

	std::string con_names[3] = {"m_0","m_Pi","mIso"};
	_con_names = con_names;
};

template <typename xdouble>
std::complex<xdouble> t_dependent_background::template_eval(const xdouble* var, const xdouble* par, const xdouble* con)const{
		xdouble m     = var[0];
		xdouble tPrime= var[1];

		xdouble b     = par[0];
		xdouble c0    = par[1];
		xdouble c1    = par[2];
		xdouble c2    = par[3];

		xdouble m0    = con[0];
		xdouble mPi   = con[1];
		xdouble mIso  = con[2];


		xdouble PSQ = 0.;
		xdouble mpor = mPi + mIso;		
		if (m > mpor){
			xdouble E = (m*m +mPi*mPi - mIso*mIso)/(2*m);
			PSQ = E*E - mPi*mPi;
		};
		return std::complex<xdouble>(pow(m-0.5,b)*exp(PSQ*(c0+c1*tPrime+c2*tPrime*tPrime)),0.);
};
// End of special definitions
#endif//BREIT_WIGNERS_WEIT_BRIGNERS
