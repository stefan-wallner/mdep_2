#ifndef MINIMIZE_MINI_MICE
#define MINIMIZE_MINI_MICE
#include<vector>
#include<complex>
#include<string>

#include"method.h"
#include"anchor_t.h"
#include"full_covariance.h"
#include"old_method.h"

typedef method METHOD; // Otherwise, the class method 'method()' and the type 'method' would have the same name --> typedef method METHOD

#include"Math/Minimizer.h"
#include"Math/Factory.h"
#include"Math/Functor.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML
double MIN_STEP_SIZE = 0.00001;

class minimize{
	public:
		minimize();
		~minimize(){delete _method;};
#ifdef USE_YAML
		minimize(std::string card);
#endif//USE_YAML
																			// In python
		double 			operator()			()				{return (*_method)(&_best_par[0]);};	
		double 			operator()			(std::vector<double>&xx)	{return (*_method)(xx);};
		double 			operator()			(const double*xx)		{return (*_method)(xx);};

		// Fitting routines
		double 			fit				();
		void 			initCouplings			(size_t nSeeds = 1		);

	// Setters and getters
		METHOD*			method()							{return _method;};
		void 			setParameter			(int i, double par		);
		void 			setParameter			(std::string name, double par	);
		void 			setParameters			(std::vector<double> pars	);
		void 			setStepSize			(int i, double step		);
		void 			setStepSize			(std::string name, double par	);
		void 			setStepSizes			(std::vector<double> steps	);
		void 			setRandRange			(double range			);
		
		double			getParameter			(size_t i)	const;			

	// Fixing and releasing
		void 			relPar(int i);
		void 			fixPar(int i);
		void 			relPar(std::string name);
		void 			fixPar(std::string name);
		void            relBranchings();
		void            fixBranchings();
		const std::vector<bool>*getReleased			()		const		{return &_released;};


	// Print routines
		std::string 		className			()		const		{return "minimize";};				// X
		void 			printStatus();													// X

	// Internal handlers
		void 			update_definitions		();
		void 			reload_par_definitions		(int mara_peter = -1);
		bool 			initialize			(std::string s1="Minuit2", std::string s2="Migrad");
		void 			setRandomCpl			();										// X
		void 			setRandomBra			();										// X
		void			findRandRange			();
		void 			finish_setUp			();

	// MULTINEST
		void			cube				(double* in)	const;

#ifdef USE_YAML	
		void 			loadFitterDefinitions(YAML::Node &waveset);
		size_t			get_method(YAML::Node &card)			const;

		YAML::Node      get_parameter_node() const;
		void            set_param_from_node(const YAML::Node& parameters);
		void            writeParamToYamlFile(const std::string& filename)const;
		void            loadParamFromYamlFile(const std::string& filename);
		void            writeResultToYamlFile(const std::string& filename)const;
#endif//USE_YAML
	protected:
	//METHOD
		METHOD*			_method;						// The method used (at the moment anchor_t)

		size_t			_method_type;
	// OWN STUFF
		std::vector<double> 	_best_par; 						// Best paramters
		double 			_randRange; 						// Range for random paramters (couplings and branchings)
		double 			_minStepSize;						// Minimal step size

	// MINIMIZER STUFF
		ROOT::Math::Minimizer* 	_min;							// ROOT Minimizer
		ROOT::Math::Functor 	_f;							// ROOT Functor object
		bool 			_init; 							// Flag for the initialization of the minimizer
		size_t 			_maxFunctionCalls;					// Miminizer definition
		size_t 			_maxIterations;						// Miminizer definition
		double 			_tolerance;						// Miminizer definition
		std::vector<double> 	_step_sizes;						// Step Size for each paramter	
		std::vector<bool> 	_released; 						// Status of each paramters

};
#endif//MINIMIZE_MINI_MICE
