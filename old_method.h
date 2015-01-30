#ifndef OLDA___METHODA
#define OLDA___METHODA
#include"waveset.h"
#include<complex>
#include<vector>
#include<string>

#include"matrix_utilities.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML

class old_method{
	public:
	// CONSTRUCTOR
		old_method();
#ifdef USE_YAML
		old_method(
								std::string 						card);
#endif//USE_YAML
	// EVALUATION METHODS
		double operator()();

		double operator()(
								std::vector<double> 					&xx);

		double operator()(
								const double						*xx);

		template<typename xdouble>
		xdouble EvalCP(
								const std::complex<xdouble>	 			*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		xdouble EvalBranch(
								const std::complex<xdouble>				*branch,
								const std::complex<xdouble>	 			*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		double EvalTbin(
								int 							tbin,
								const std::complex<double> 				*cpl,
								const double	 					*par,
								const double						*iso_par)			const;
#ifdef ADOL_ON
		adouble EvalTbin(
								int 							tbin,
								const std::complex<adouble> 				*cpl,
								const adouble	 					*par,
								const adouble						*iso_par)			const;
#endif//ADOL_ON


		template<typename xdouble>
		xdouble EvalBin(
								int 							tbin,
								int 							bin,
								const std::complex<xdouble> 				*cpl,
								const xdouble	 					*par,
								std::vector<std::vector<std::complex<xdouble> > > 	&iso_eval)			const;

		template<typename xdouble>
		std::vector<xdouble> delta(
								int 							tbin,
								int 							bin,
								double 							mass,
								const std::complex<xdouble>	 			*cpl,
								const xdouble	 					*par,
								std::vector<std::vector<std::complex<xdouble> > > 	&iso_eval)			const;

	// DERIVATIVES
#ifdef ADOL_ON
		std::vector<double> 				Diff(std::vector<double> &xx)								const;
		std::vector<double> 				Diff(const double* xx)									const;
#endif//ADOL_ON
	// PARAMETER HANDLING
		void 					setParameter(size_t i, double par);
		bool 					setParameter(std::string name, double par);
		void 					setParameters(std::vector<double> pars);
		int 					getParNumber(std::string name)		const;
		void 					setParLimits(int i, double upper, double lower);
		void					init_lower_limits(int n=-1);
		void					init_upper_limits(int n=-1);
		std::vector<std::complex<double> > 	get_branchings(const std::vector<std::complex<double> > &cpl,const std::vector<double> &par,const std::vector<double> &iso_par) const;
		std::vector<std::complex<double> >	getUnbranchedCouplings(const std::vector<std::complex<double> > &cpl,const std::vector<std::complex<double> > &bra) const;
		std::vector<std::complex<double> >	getAllCouplings(int tbin,const std::vector<std::complex<double> > &cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const;
		void 					branchCouplingsToOne();



	// DATA HANDLING
		bool 					set_data(int tbin, int bin, std::vector<double> data);
		bool 					set_coma(int tbin, int bin, std::vector<double> coma);
		void 					loadData(int tbin, const char* dataFile);
		void 					loadComa(int tbin, const char* comaFile);
		void 					nullify();
		void 					conjugate();
#ifdef USE_YAML
	// YAML SETTER
		bool					loadDataComa(YAML::Node &waveset);
		bool					loadParameterValues(YAML::Node &waveset, YAML::Node &param);
#endif//USE_YAML
	// OTHER SETTERS & GETTERS
		waveset* 				Waveset			()		{return &_waveset;};
		size_t					nTot			()const		{return _nTot;};
		size_t					nPar			()const		{return _nPar;};
		size_t					nCpl			()const		{return _nCpl;};
		size_t					nBra			()const		{return _nBra;};
		size_t 					nIso			()const		{return _nIso;};
		size_t					nBrCpl			()const		{return _nBrCpl;};
		const std::vector<double>		parameters		()const;
		const std::vector<double>*		upper_parameter_limits	()const		{return &_upper_parameter_limits;};
		const std::vector<double>*		lower_parameter_limits	()const		{return &_lower_parameter_limits;};
		const std::vector<double>* 		get_data(int tbin, int bin)const	{return &_data[tbin][bin];};
		const std::vector<double>*		get_coma(int tbin, int bin)const	{return &_coma[tbin][bin];};

		const std::vector<std::string>*		parNames		()const		{return &_parNames;};

		std::vector<double>			getParameters		()const;
		double					getParameter(size_t i)	const;

		int getNtot();
		int getNcpl();

		std::vector<std::vector<std::complex<double> > > full_to_br_cpl(std::vector<std::complex<double> > &cpl);

	// OTHER METHODS
		std::string 				className		()const		{return "old_method";};
		void 					printStatus()		const;
		void 					setTbinning(std::vector<std::vector<double> > binning);
		void 					update_n_cpls();
		void 					update_min_max_bin();
		void 					update_definitions();
		void					update_is_active();

	// PLOTTING
		void					write_plots(std::string filename, int tbin) const;
		void					write_plots(std::string filename, int tbin,const std::vector<double> &paramters) const;
		void					write_plots(std::string filename, int tbin,const std::vector<std::complex<double> >&cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const;

	protected:
		// WAVESET
		waveset									_waveset;		// The waveset used

		// PARAMETER NUMBERS

		size_t 									_nTot; 			// Total number of parameters
		size_t 									_nPar; 			// Number of shape parameters
		size_t 									_nCpl; 			// Number of couplings (total, all t' bins summed)
		size_t 									_nBra; 			// Number of branchings
		size_t 									_nIso;			// Number of isobar parameters
		size_t 									_nBrCpl;		// Number of couplings with branchings in the anchor wave

		// PARAMETERS AND DATA
		std::vector<double>							_parameters; 		// Acutal paramters (2*_nCpl,_nPar,2*_nBra,_nIso) - these are 'all' parameters!!!
		std::vector<double> 							_upper_parameter_limits;// Paramters limits
		std::vector<double> 							_lower_parameter_limits;
		std::vector<std::string> 						_parNames; 		// Name of each parameter
		std::vector<std::vector<std::vector<double> > > 			_data; 			// Data
		std::vector<std::vector<std::vector<double> > > 			_coma; 			// Just the errors in this case
		std::vector<std::vector<bool> >						_is_active;		// Flag, which point is actually active

		// OTHER MEMBERS
		size_t 									_nOut; 			// Print output after _nOut iterations
		size_t 									_count;			// Count of calls
};


#endif//OLDA___METHODA
