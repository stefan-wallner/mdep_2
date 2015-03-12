#ifndef METHOD_METHOD_METHOD
#define METHOD_METHOD_METHOD
#include"waveset.h"
#include<string>
#include<stdexcept>


class method{
	public:
		method():_waveset(),_nOut(1000),_count(0){};
		method(std::string card):_waveset(card),_nOut(1000),_count(0){};

		double				operator()						()											{return mainEval(&parameters()[0]);};
		double 				operator()						(std::vector<double> &xx)								{return mainEval(&xx[0]);};
		double 				operator()						(const double *xx)									{return mainEval( xx);};

		virtual double 			mainEval						(const double *xx)									{std::cout<<"method.h: mainEval(...) not overwritten"<<std::endl;throw;return 0.;};

	// OTHER SETTERS & GETTERS
		waveset* 				Waveset						()											{return &_waveset;};
		size_t					nTot						()										const	{return _nTot;};
		size_t					nPar						()										const	{return _nPar;};
		size_t					nCpl						()										const	{return _nCpl;};
		size_t					nBra						()										const	{return _nBra;};
		size_t 					nIso						()										const	{return _nIso;};
		size_t					nBrCpl						()										const	{return _nBrCpl;};
		const std::vector<double>*		upper_parameter_limits				()										const	{return &_upper_parameter_limits;};
		const std::vector<double>*		lower_parameter_limits				()										const	{return &_lower_parameter_limits;};
		const std::vector<double>* 		get_data					(int tbin, int bin)								const	{return &_data[tbin][bin];};
		const std::vector<std::string>*		parNames					()										const	{return &_parNames;};

		void 					setParameter					(size_t i, double par);
		double					getParameter					(size_t i)									const;
		void 					setParameters					(std::vector<double> pars);
		const std::vector<double>		parameters					()const;
		bool 					setParameter					(std::string name, double par);
		int 					getParNumber					(std::string name)								const;
		void 					setParLimits					(int i, double upper, double lower);
		void					init_lower_limits				(int n=-1);
		void					init_upper_limits				(int n=-1);
		std::vector<std::complex<double> >	getUnbranchedCouplings				(const std::vector<std::complex<double> > &cpl,const std::vector<std::complex<double> > &bra) 	const;
		void 					update_min_max_bin				();

		virtual std::string			className					()										const	{return "undefined_method";};

		virtual void 				update_definitions				()											{std::cout<<"method.h: update_definitions() not overwritten"<<std::endl;throw;};
		virtual void 				update_n_cpls					()											{std::cout<<"method.h: update_n_cpls() not overwritten"<<std::endl;throw;};

		virtual void 				printStatus					()										const	{std::cout<<"method.h: printStatus() not overwritten"<<std::endl;throw;};

		virtual void 				conjugate					()											{std::cout<<"method.h: conjugate() not overwritten"<<std::endl;throw;};

		virtual void 				loadData					(int tbin, const char* dataFile)							{std::cout<<"method.h: loadData(...) not overwritten"<<std::endl;throw;};

		virtual void 				loadComa					(int tbin, const char* comaFile)							{std::cout<<"method.h: loadComa(...) not overwritten"<<std::endl;throw;};

// This is a special method. The base classes just return 'false' autoCpl classes return 'true' and do something
		virtual bool 				setUseBranch					(bool in)										{return false;};

#ifdef USE_YAML
	// YAML SETTER
		bool					loadDataComa					(YAML::Node &waveset);
		bool					loadParameterValues				(YAML::Node &waveset, YAML::Node &param);

		YAML::Node      get_parameter_node() const;
		void            set_param_from_node(const YAML::Node& parameters);
		void            writeParamToFile(const std::string& filename)const;
		void            loadParamFromFile(const std::string& filename);
#endif//USE_YAML

	// PLOTTING
		void					write_plots					(std::string filename, int tbin) 						const;
		void					write_plots					(std::string filename, int tbin,const std::vector<double> &paramters) 		const;
		virtual void				write_plots(std::string filename, int tbin,const std::vector<std::complex<double> >&cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso) const {std::cout<<"method.h: write_plots(...) not overwritten"<<std::endl;throw;};
		virtual std::vector<std::complex<double> > getAllCouplings(int tbin,const std::vector<std::complex<double> > &cpl,const std::vector<double> &par,const std::vector<std::complex<double> > &bra,const std::vector<double> &iso)const{std::cout<<"method.h: getAllCouplings(...) not overwritten"<<std::endl;throw;};

		std::vector<std::complex<double> > amplitudes					(double m, int tbin) const;
		std::vector<std::complex<double> > amplitudes					(double m, int tbin, const std::vector<double> &parameters) const;

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

		// OTHER MEMBERS
		size_t 									_nOut; 			// Print output after _nOut iterations
		size_t 									_count;			// Count of calls

};
#endif//METHOD_METHOD_METHOD
