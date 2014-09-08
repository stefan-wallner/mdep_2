#ifndef ANCHOR_ROHCNA
#define ANCHOR_ROHCNA
#include"waveset.h"
#include<complex>
#include<vector>
#include<string>
#include"AB.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML

class anchor_t : public waveset {
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
								std::vector<std::complex<xdouble> > 			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);

		template<typename xdouble>
		xdouble EvalBranch(
								std::vector<std::complex<xdouble> >			&branch,
								std::vector<std::complex<xdouble> > 			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);

		template<typename xdouble>
		xdouble EvalAutoCpl(
								std::vector<std::complex<xdouble> > 			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);

		template<typename xdouble>
		xdouble EvalAutoCplBranch(
								std::vector<std::complex<xdouble> >			&bra,
								std::vector<std::complex<xdouble> >			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);

		template<typename xdouble>
		xdouble EvalAutoCplTbin(
								int 							tbin,
								std::vector<std::complex<xdouble> > 			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);

		template<typename xdouble>
		xdouble EvalTbin(
								int 							tbin,
								std::vector<std::complex<xdouble> > 			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);

		template<typename xdouble>
		xdouble EvalBin(
								int 							tbin,
								int 							bin,
								std::vector<std::complex<xdouble> > 			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<std::vector<std::complex<xdouble> > > 	&iso_eval);

		template<typename xdouble>
		std::vector<xdouble> delta(
								int 							tbin,
								int 							bin,
								double 							mass,
								std::vector<std::complex<xdouble> > 			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<std::vector<std::complex<xdouble> > > 	&iso_eval);

	// AUTO CPL METHODS
		template<typename xdouble>
		AandB<xdouble> get_AB(
								int 							tbin,
								std::vector<std::complex<xdouble> > 			&anchor_cpl,
								std::vector<xdouble> 					&par);

		template<typename xdouble>
		AandB<xdouble> get_AB_iso( // For de-isobarred
								int 							tbin,
								std::vector<std::complex<xdouble> > 			&anchor_cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);

		template<typename xdouble>
		std::vector<std::complex<xdouble> > getMinimumCpl(
								int 							tbin,
								std::vector<std::complex<xdouble> > 			&anchor_cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);

		template<typename xdouble>
		std::vector<std::complex<xdouble> > getMinimumCplBra(
								int 							tbin,
								std::vector<std::complex<xdouble> > 			&branch,
								std::vector<std::complex<xdouble> > 			&anchor_cpl,
								std::vector<xdouble> 					&par,
								std::vector<xdouble> 					&iso_par);
	// DERIVATIVES
#ifdef ADOL_ON
		std::vector<double> 				Diff(std::vector<double> &xx);
		std::vector<double> 				Diff(const double* xx);
#endif//ADOL_ON
	// PARAMETER HANDLING
		void 					setParameter(int i, double par);
		bool 					setParameter(std::string name, double par);
		void 					setParameters(std::vector<double> pars);
		std::vector<double> 			getParameters();
		std::vector<std::string>		getParNames();
		std::string 				getParName(int i);
		int 					getParNumber(std::string name);
		void 					setParLimits(int i, double upper, double lower);
		std::vector<std::complex<double> > 	get_branchings(std::vector<std::complex<double> > &cpl,std::vector<double> &par, std::vector<double> &iso_par);
		std::vector<std::complex<double> >	getUnbranchedCouplings(std::vector<std::complex<double> > &cpl, std::vector<std::complex<double> > &bra);
		std::vector<std::complex<double> >	getAllCouplings(int tbin,std::vector<std::complex<double> > &cpl, std::vector<double> &par, std::vector<std::complex<double> > &bra, std::vector<double> &iso);
		void 					branchCouplingsToOne();



	// DATA HANDLING
		bool 					set_data(int tbin, int bin, std::vector<double> data);
		bool 					set_coma(int tbin, int bin, std::vector<std::vector<double> > coma);
		std::vector<double> 			get_data(int tbin, int bin);
		std::vector<std::vector<double> > 	get_coma(int tbin, int bin);
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
		int getNtotAnc();
		int getNanc();

	// OTHER METHODS
		std::string 				className();
		void 					printStatus();
		void 					set_is_ampl(bool is_ampl);
		void 					setTbinning(std::vector<double> binning);
		void 					update_n_cpls();
		void 					update_min_max_bin();
		void 					update_definitions();

	// PLOTTING
		void					write_plots(std::string filename, int tbin);
		void					write_plots(std::string filename, int tbin, std::vector<double> &paramters);
		void					write_plots(std::string filename, int tbin,std::vector<std::complex<double> >&cpl,std::vector<double> &par, std::vector<std::complex<double> > &bra, std::vector<double> &iso);

	protected:
		// PARAMETER NUMBERS
		int 									_nTot; 			// Total number of parameters
		int 									_nPar; 			// Number of shape parameters
		int 									_nCpl; 			// Number of couplings (total, all t' bins summed)
		int 									_nBra; 			// Number of branchings
		int 									_nIso;			// Number of isobar parameters
		int 									_nBrCplAnc;		// Number of couplings with branchings in the anchor wave

		// PARAMETERS AND DATA
		std::vector<double>							_parameters; 		// Acutal paramters (2*_nCpl,_nPar,2*_nBra) - these are 'all' parameters!!!
		std::vector<double> 							_upper_parameter_limits;// Paramters limits
		std::vector<double> 							_lower_parameter_limits;
		std::vector<std::string> 						_parNames; 		// Name of each parameter
		std::vector<std::vector<std::vector<double> > > 			_data; 			// Data
		std::vector<std::vector<std::vector<std::vector<double> > > > 		_coma; 			// Covariance matrix

		// OTHER MEMBERS
		bool 									_is_ampl; 		// Enable Amplitunde fitting, needs to be tested.
		bool 									_useBranch; 		// Switches the usage of branchings on/off (only in the operator() method
		int 									_nOut; 			// Print output after _nOut iterations
		int 									_count;			// Count of calls
};


#endif//ANCHOR_ROHCNA
