#ifndef MINIMIZE_MINI_MICE
#define MINIMIZE_MINI_MICE
#include<vector>
#include<complex>
#include<string>
#include "anchor_t.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

class minimize : public anchor_t{
	public:
		minimize();

		double operator()(const double* xx);
		double operator()();

		std::vector<std::string> getParNames();
		std::vector<bool> getReleased();
		std::vector<double> getParameters();
		void setParameter(std::string name, double par);
		void setStepSize(std::string name, double par);
		void fixPar(std::string name);
		void relPar(std::string name);
		int getParNumber(std::string name);
		void setParameter(int i, double par);
		void setParameters(std::vector<double> pars);
		void setParLimits(int i, double upper, double lower);
		void setStepSize(int i, double step);
		void setStepSizes(std::vector<double> steps);
		std::string getParName(int i);
		void relPar(int i);
		void fixPar(int i);

		void reload_par_definitions(int mara_peter = -1);
		bool initialize(std::string s1="Minuit2", std::string s2="Migrad");	
		void update_definitions();
		void update_parameters();
		void setRandomCpl();
		void setRandomBra();
		void setRandRange(double range);
		void initCouplings();
		void printStatus();

		double fit();

		void writePlots(int tbin,std::string filename);
		void branchCouplingsToOne();
	protected:
		bool _init; // true, if the minimizer is already initialized
		bool _useBranch; // Switches the usage of branchings on/off (only in the operator() method

		int _nOut; // Print output after _nOut iterations
		int _nTot; // Total number of parameters
		int _nPar; // Number of shape parameters
		int _nCpl; // Number of couplings (total, all t' bins summed)
		int _nBra; // Number of branchings
		int _count;// Count of calls
		int _maxFunctionCalls;
		int _maxIterations;

		double _tolerance;
		double _randRange; // Range for random paramters (couplings and branchings)
		std::vector<double> _best_par; // Best paramters
		ROOT::Math::Minimizer* _min;
		ROOT::Math::Functor _f;
		double _minStepSize; // Minimal step size

		std::vector<double> _parameters; // Acutal paramters (2*_nCpl,_nPar,2*_nBra) - these are 'all' parameters!!! In the lower classes, parameters were only the shape parameters
		std::vector<double> _upper_parameter_limits; // Paramters limits
		std::vector<double> _lower_parameter_limits;
		std::vector<double> _step_sizes;
		std::vector<bool> _released; // status of each paramters
		std::vector<std::string> _parNames; // Name of each parameter
};
#endif//MINIMIZE_MINI_MICE
