#ifndef ANCHOR_ROHCNA
#define ANCHOR_ROHCNA
#include"waveset.h"
#include<complex>
#include<vector>
#include<string>

#include"matrix_utilities.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML

class anchor_t{
	public:
	// CONSTRUCTOR
		anchor_t();
#ifdef USE_YAML
		anchor_t(
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

		template<typename xdouble>
		xdouble EvalAutoCpl(
								const std::complex<xdouble>	 			*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		xdouble EvalAutoCplBranch(
								const std::complex<xdouble>				*bra,
								const std::complex<xdouble>				*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		xdouble EvalAutoCplTbin(
								int 							tbin,
								const std::complex<xdouble> 				*cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		xdouble EvalTbin(
								int 							tbin,
								const std::complex<xdouble> 				*cpl,
								const xdouble	 					*par,
								const xdouble						*iso_par)			const;

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

	// AUTO CPL METHODS
		template<typename xdouble>
		AandB<xdouble> get_AB( // For de-isobarred
								int 							tbin,
								const std::complex<xdouble>	 			*anchor_cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		std::vector<std::complex<xdouble> > getMinimumCpl(
								int 							tbin,
								const std::complex<xdouble> 				*anchor_cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;

		template<typename xdouble>
		std::vector<std::complex<xdouble> > getMinimumCplBra(
								int 							tbin,
								const std::complex<xdouble>	 			*branch,
								const std::complex<xdouble> 				*anchor_cpl,
								const xdouble	 					*par,
								const xdouble	 					*iso_par)			const;
	// DERIVATIVES
#ifdef ADOL_ON
		std::vector<double> 				Diff(std::vector<double> &xx)								const;
		std::vector<double> 				Diff(const double* xx)									const;
#endif//ADOL_ON
	// PARAMETER HANDLING
		void 					setParameter(int i, double par);
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
		bool 					set_coma(int tbin, int bin, std::vector<std::vector<double> > coma);
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
		bool					useBranch		()const		{return _useBranch;};
		int					nTot			()const		{return _nTot;};
		int					nPar			()const		{return _nPar;};
		int					nCpl			()const		{return _nCpl;};
		int					nBra			()const		{return _nBra;};
		int 					nIso			()const		{return _nIso;};
		int					nBrCplAnc		()const		{return _nBrCplAnc;};
		const std::vector<double>		parameters		()const;
		const std::vector<double>*		upper_parameter_limits	()const		{return &_upper_parameter_limits;};
		const std::vector<double>*		lower_parameter_limits	()const		{return &_lower_parameter_limits;};
		const std::vector<double>* 		get_data(int tbin, int bin)const	{return &_data[tbin][bin];};
		const std::vector<std::vector<double> >*get_coma(int tbin, int bin)const	{return &_coma[tbin][bin];};

		const std::vector<std::string>*		parNames		()const		{return &_parNames;};

		std::vector<double>			getParameters		()const;
		double					getParameter(size_t i)	const;

		int getNtotAnc();
		int getNanc();
		void setUseBranch(bool in);


	// OTHER METHODS
		std::string 				className		()const		{return "anchor_t";};
		void 					printStatus()		const;
		void 					set_is_ampl(bool is_ampl);
		void 					setTbinning(std::vector<double> binning);
		void 					update_n_cpls();
		void 					update_min_max_bin();
		void 					update_definitions();
#ifdef STORE_ACTIVE
		void					update_is_active();
#endif//STORE_ACTIVE

	// PLOTTING
		void					write_plots(std::string filename, int tbin) const;
		void					write_plots(std::string filename, int tbin,const std::vector<double> &paramters) const;
		void					write_plots(std::string filename, int tbin,const std::vector<std::complex<double> >&cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const;
	protected:
		// WAVESET
		waveset									_waveset;		// The waveset used

		// PARAMETER NUMBERS

		int 									_nTot; 			// Total number of parameters
		int 									_nPar; 			// Number of shape parameters
		int 									_nCpl; 			// Number of couplings (total, all t' bins summed)
		int 									_nBra; 			// Number of branchings
		int 									_nIso;			// Number of isobar parameters
		int 									_nBrCplAnc;		// Number of couplings with branchings in the anchor wave

		// PARAMETERS AND DATA
		std::vector<double>							_parameters; 		// Acutal paramters (2*_nCpl,_nPar,2*_nBra,_nIso) - these are 'all' parameters!!!
		std::vector<double> 							_upper_parameter_limits;// Paramters limits
		std::vector<double> 							_lower_parameter_limits;
		std::vector<std::string> 						_parNames; 		// Name of each parameter
		std::vector<std::vector<std::vector<double> > > 			_data; 			// Data
		std::vector<std::vector<std::vector<std::vector<double> > > > 		_coma; 			// Covariance matrix
#ifdef STORE_ACTIVE
		std::vector<std::vector<std::vector<bool> > >				_is_active;		// Flag, which point is actually active
#endif//STORE_ACTIVE

		// OTHER MEMBERS
		bool 									_is_ampl; 		// Enable Amplitunde fitting, needs to be tested.
		bool 									_useBranch; 		// Switches the usage of branchings on/off (only in the operator() method
		int 									_nOut; 			// Print output after _nOut iterations
		int 									_count;			// Count of calls
};


#endif//ANCHOR_ROHCNA
