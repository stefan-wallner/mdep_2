#ifndef MINIMIZE_MINI_MICE
#define MINIMIZE_MINI_MICE
#include<vector>
#include<complex>
#include<string>
#include"anchor_t.h"
#include"Math/Minimizer.h"
#include"Math/Factory.h"
#include"Math/Functor.h"

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML
double MIN_STEP_SIZE = 0.00001;

class minimize : public anchor_t{
	public:
		minimize();

#ifdef USE_YAML
		minimize(std::string card);
#endif//USE_YAML
	// Fitting routines
		double 			fit();
		void 			initCouplings();

	// Setters and getters
		void 			setParameter(int i, double par);
		void 			setParameter(std::string name, double par);
		void 			setParameters(std::vector<double> pars);
		void 			setStepSize(int i, double step);
		void 			setStepSize(std::string name, double par);
		void 			setStepSizes(std::vector<double> steps);
		void 			setRandRange(double range);

	// Fixing and releasing
		void 			relPar(int i);
		void 			fixPar(int i);
		void 			relPar(std::string name);
		void 			fixPar(std::string name);
		std::vector<bool> 	getReleased();


	// Print routines
		std::string 		className();
		void 			printStatus();

	// Internal handlers
		void 			update_definitions();
		void 			reload_par_definitions(int mara_peter = -1);
		bool 			initialize(std::string s1="Minuit2", std::string s2="Migrad");
		void 			setRandomCpl();
		void 			setRandomBra();
		void 			finish_setUp();

#ifdef USE_YAML	
		void 			loadFitterDefinitions(YAML::Node &waveset);
#endif//USE_YAML
	protected:
	// OWN STUFF
		std::vector<double> 	_best_par; 						// Best paramters
		double 			_randRange; 						// Range for random paramters (couplings and branchings)
		double 			_minStepSize;						// Minimal step size

	// MINIMIZER STUFF
		ROOT::Math::Minimizer* 	_min;							// ROOT Minimizer
		ROOT::Math::Functor 	_f;							// ROOT Functor object
		bool 			_init; 							// Flag for the initialization of the minimizer
		int 			_maxFunctionCalls;					// Miminizer definition
		int 			_maxIterations;						// Miminizer definition
		double 			_tolerance;						// Miminizer definition
		std::vector<double> 	_step_sizes;						// Step Size for each paramter	
		std::vector<bool> 	_released; 						// Status of each paramters


};
#endif//MINIMIZE_MINI_MICE
