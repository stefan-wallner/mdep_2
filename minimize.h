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
		minimize(std::string card, std::string waves, std::string parametrizations);
#endif//USE_YAML
		void setStepSize(std::string name, double par);

		void setParameter(int i, double par);
		void setParameter(std::string name, double par);
		void setParameters(std::vector<double> pars);
		std::string className();
		std::vector<bool> getReleased();
		void fixPar(std::string name);
		void relPar(std::string name);
		void setStepSize(int i, double step);
		void setStepSizes(std::vector<double> steps);
		void relPar(int i);
		void fixPar(int i);
		void reload_par_definitions(int mara_peter = -1);
		bool initialize(std::string s1="Minuit2", std::string s2="Migrad");
		void initCouplings();
		void printStatus();
		double fit();
		void update_definitions();
		void setRandomCpl();
		void setRandomBra();
		void setRandRange(double range);
		void finish_setUp();
#ifdef USE_YAML
		void loadFitterDefinitions(YAML::Node &waveset);
#endif//USE_YAML
	protected:
		bool _init; // true, if the minimizer is already initialized

		int _maxFunctionCalls;
		int _maxIterations;

		double _tolerance;

		std::vector<double> _best_par; // Best paramters
		ROOT::Math::Minimizer* _min;
		ROOT::Math::Functor _f;
		double _minStepSize; // Minimal step size
		double _randRange; 		// Range for random paramters (couplings and branchings)
		std::vector<double> _step_sizes;
		std::vector<bool> _released; // status of each paramters

};
#endif//MINIMIZE_MINI_MICE
