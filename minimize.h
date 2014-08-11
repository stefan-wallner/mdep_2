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



	protected:
		int _n_self;

		bool _init;
		bool _useBranch;
		int _nOut;

		int _nTot;
		int _nPar;
		int _nCpl;
		int _nBra;
		int _count;

		int _maxFunctionCalls;
		int _maxIterations;
		double _tolerance;
		double _randRange;

		std::vector<double> _best_par;

		ROOT::Math::Minimizer* _min;
		ROOT::Math::Functor _f;

		double _minStepSize;
		std::vector<double> _parameters;
		std::vector<double> _upper_parameter_limits;
		std::vector<double> _lower_parameter_limits;

		std::vector<double> _step_sizes;
		std::vector<bool> _released;
		std::vector<std::string> _parNames;

};
#endif//MINIMIZE_MINI_MICE
