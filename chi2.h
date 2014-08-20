#ifndef CHI2_LOLOLO
#define CHI2_LOLOLO
#include<vector>
#include<complex>
#include<fstream>
#include<iostream>

class chi2{
	public:
		// CONSTRUCTORS
		chi2();
		// EVALUATE THE FUNCTION

		template<typename xdouble>
		std::vector<std::complex<xdouble> > amps(double m,std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par);

		template<typename xdouble>
		std::vector<std::complex<xdouble> > funcs(double m,std::vector<xdouble> &par);
		std::vector<double> phase_space(double m);
		// MANAGE FUNCTIONS
		std::string className();
		void add_wave();
		void add_func(int i);
		void add_func_to_wave(int wave, int func);

		// HELPERS
		bool checkConsistency();

		std::vector<int> get_wave_functions(int wave);
		std::vector<int> get_wave_pars(int wave);
		std::vector<int> get_wave_const(int wave);
		std::vector<int> get_function_pars(int func);
		std::vector<int> get_function_const(int func);
		std::vector<int> get_function_waves(int func);

		std::vector<int> get_nPars();
		std::vector<int> get_nConst();

		void setWaveName(int i, std::string name);
		void setFunctionName(int i, std::string name);
		void setParameterName(int i, std::string name);
		void setConstantName(int i, std::string name);
		std::string getWaveName(int i);
		std::string getFunctionName(int i);
		std::string getParameterName(int i);
		std::string getConstantName(int i);

		void setWaveLimits(int i, double lower, double upper);
		void setWaveSpin(int i, int L);
		void setWavePhaseSpace(int i, int ps);
		void setConst(int i,double con);
		void setGlobalPhaseSpace(int i);
		void updateFuncLims();
		void updateFuncSpin();
		void updateNftw();

		void printStatus();
		void printParameters();

		// METHODS FOR INTERFACING WITH THE 'OLD' PARAMETERS [RE,IM,PAR,...,RE,IM,PAR...]
		void setInterface();
		std::vector<std::complex<double> > amps_class(double m, std::vector<double> &par);
		std::vector<std::complex<double> > cpls(std::vector<double> &par);
		std::vector<double> pars(std::vector<double> &par);

		void open_output(std::string filename ="chi2.dat");
		void close_output();

	protected:

		int _nWaves; 					// Number of waves
		int _nFuncs; 					// Number of functions

		int _maxNpars;					// maximum Number of parameters a BW function has
		int _nPar; 					// Number of parameters
		int _nFtw; 					// Number of function-wave couplings (should be _funcs_to_waves.size())
		std::vector<int> _funcs; 			// Defined functions (the int gives the number of the parametrization in 'breitWigners.h')
		std::vector<int> _borders_waves; 		// Tells, which entries in funcs to waves belong to which wave
		std::vector<int> _funcs_to_waves; 		// Tells, which function is to be coupled to the wave (the blocks for each wave are encoded in '_borders_waves')
		std::vector<int> _borders_par; 			// Tells, which paramters belong to which function
		std::vector<int> _borders_const;		// Tells, which constants belong to which function
		std::vector<int> _L; 				// Gives the Spin of each wave
		std::vector<int> _L_func; 			// Gives the spin of for each function
		std::vector<double> _upperLims; 		// Gives     limits	each
		std::vector<double> _lowerLims; 		// 	 the	    for	     wave
		std::vector<double> _funcUpperLims; 		// Mass limits for each function
		std::vector<double> _funcLowerLims; 		// ''	''	''	''
		std::vector<double> _const; 			// Values for the constants in the functions

		std::vector<int> _interface;			// Some inteface to be able to read parameter stes from the old version (not really needed anymore)

		// Phase space definitions
		int _globalPs; 					// Assume only ONE global PS and   //Also assume ps that does not depend on the isobar mass
		std::vector<int> _wavePs; 			// One PS for each wave // Or the spin

		// Names
		std::vector<std::string> _waveNames;
		std::vector<std::string> _funcNames;
		std::vector<std::string> _parNames;
		std::vector<std::string> _constNames;

		// Analysis stuff
		bool _write_out;				// Flag to switch on the text_output
		std::ofstream *_outStream;			// Stream for the text output



};
#endif//CHI2_LOLOLO
