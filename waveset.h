#ifndef WAVESET_SETS_WAVES
#define WAVESET_SETS_WAVES
#include<vector>
#include<complex>
#include<fstream>
#include<iostream>


class waveset {
	public:
	// CONSTRUCTOR
		waveset();

	// FUNCTIONS AND AMPLITUDES

		template<typename xdouble>
		std::vector<std::complex<xdouble> > amps(	
								double 							m,
								std::vector<std::complex<xdouble> > 			&cpl,
								std::vector<xdouble> 					&par,
								std::vector<std::vector<std::complex<xdouble> > > 	&funcEvals2pi);

		template<typename xdouble>
		std::vector<std::complex<xdouble> > funcs(	
								double 							m,
								std::vector<xdouble> 					&par);

		template<typename xdouble>
		std::vector<std::vector<std::complex<xdouble> > > iso_funcs(
								std::vector<xdouble> 					&par);

		std::vector<double> phase_space(
								double 							m);

	// SET UP WAVESET
		// // WAVES AND FUNCTIONS
		void add_wave();
		void add_func(int i);
		void add_iso(int i);

		// // SET AMPLITUDE DEFINITIONS
		void add_func_to_wave(int wave, int func);
		void add_funcs_to_wave(int wave, int func, int func2);

		// // SETTINGS FOR DIFFERENT WAVES
		void setWaveLimits(int i, double lower, double upper);
		void setWaveSpin(int i, int L);
		void setWaveIsobarSpin(int wave, int L);
		void setGlobalPhaseSpace(int i);
		void setWavePhaseSpace(int i, int ps);
		void setWaveIsobarBinning(int wave, int binning);
		void setConst(int i,double con);
		void set_iso_const(int con, double value);
		void add_isobar_binning(std::vector<double> binning);

		// // NAMES
		void setWaveName(int i, std::string name);
		void setFunctionName(int i, std::string name);
		void setParameterName(int i, std::string name);
		void setConstantName(int i, std::string name);
		void setIsobarName(int i,std::string name);
		void setIsoParName(int i,std::string name);
		void setIsoConstName(int i,std::string name);

	// GETTERS
		// // SIMPLE OVER ALL NUMBERS
		int getNpoints();
		int getNftw();		

		// // PROPERTIES OF THE WAVES
		std::string getWaveName(int i);
		std::vector<int> get_wave_functions(int wave);
		std::vector<int> get_wave_pars(int wave);
		std::vector<int> get_wave_const(int wave);
		std::vector<int> get_wave_isobars(int wave);
		std::vector<int> get_wave_iso_pars(int wave);
		std::vector<int> get_wave_iso_const(int wave);

		// // PROPERTIES OF THE FUNCTIONS
		std::string getFunctionName(int i);
		std::vector<int> get_nPars();
		std::vector<int> get_nConst();
		std::vector<int> get_function_pars(int func);
		std::vector<int> get_function_const(int func);
		std::vector<int> get_function_waves(int func);

		// // PROPERTIES OF THE ISOBARS
		std::string getIsobarName(int i);
		std::vector<int> get_nParsIso();
		std::vector<int> get_nConstIso();
		std::vector<int> get_isobar_pars(int func);
		std::vector<int> get_isobar_const(int func);
		std::vector<int> get_isobar_waves(int func);

		// // PROPERTIES OF THE PARAMETERS
		std::string getParameterName(int i);
		std::string getConstantName(int i);
		std::string getIsoParName(int i);
		std::string getIsoConstName(int i);

	// UPDATERS
		void updateNftw();
		void updateNpoints();
		void updateFuncLims();
		void updateFuncSpin();
		void updateIsobar();

	// INFO FUNCTIONS
		std::string className();
		bool checkConsistency();
		void printStatus();
		void printParameters();
		void open_output(std::string filename ="chi2log.dat");
		void close_output();

	protected:
	// WAVES
		int 					_nWaves; 		// Number of waves
		std::vector<std::string> 		_waveNames;
		int					_nPoints; 		// Number of amplitudes in the end (+1 for each wave, +1 for each isobar step)
		int 					_nFtw; 			// Number of function-wave couplings (should be _funcs_to_waves.size())
		std::vector<double>			_upperLims; 		// Gives     limits	each
		std::vector<double> 			_lowerLims; 		// 	 the	    for	     wave
		std::vector<int> 			_L; 			// Gives the Spin of each wave
		std::vector<int> 			_L_iso; 		// Gives the spin of the isobar for the single waves
		std::vector<int> 			_wave_binning_pts; 	// Number of bins for the isobars
		std::vector<int> 			_point_to_wave;		// Tells, which points belong to which wave
		std::vector<int> 			_wave_n_binning; 	// Tells, which isobar binning for the single isobars

	// FUNCTIONS
		std::vector<int> 			_funcs_to_waves; 	// Tells, which function is to be coupled to the wave (the blocks for each wave are encoded in '_borders_waves')
		int 					_nPar; 			// Number of parameters
		int 					_nFuncs; 		// Number of functions
		std::vector<std::string> 		_funcNames;
		std::vector<int> 			_funcs; 		// Defined functions (the int gives the number of the parametrization in 'breitWigners.h')
		std::vector<int> 			_borders_waves; 	// Tells, which entries in funcs to waves belong to which wave
		int 					_maxNpars;		// maximum Number of parameters a BW function has
		std::vector<double> 			_funcUpperLims; 	// Mass limits for each function
		std::vector<double> 			_funcLowerLims; 	// ''	''	''	''
		std::vector<int> 			_L_func; 		// Gives the spin of for each function

	// ISOBARS
		std::vector<int> 			_iso_to_waves; 		// Couples isobars to waves (-1 -> no isobar paramtrization -> Standard mdep fit)
		int 					_nIso; 			// Number of isobars
		std::vector<std::string> 		_iso_funcNames;
		std::vector<int> 			_isos; 			// Isobar paramterizations
		int 					_maxNparsIso; 		// Maximum number of parameters for one isobar
		int 					_maxBinIso; 		// Maximum number of bins for one isobar
		std::vector<int> 			_iso_n_binning; 	// Tells, which isobar binning to use for the waves
		std::vector<std::vector<double> >	_iso_binnings;		// Gives different isobar binnings
		std::vector<int> 			_point_borders_wave;	// Tells, where the isobar points for a single wave start
		std::vector<int> 			_L_iso_func; 		// Gives the spin of the the single isobars
		std::vector<int> 			_iso_binning_pts; 	// Number of isobar bins for each wave (no isobar -> -1)

	// PARAMETERS & CONSTANTS
		std::vector<std::string> 		_parNames;
		std::vector<int> 			_borders_par; 		// Tells, which paramters belong to which function
		std::vector<std::string>		_constNames;
		std::vector<int>			_borders_const;		// Tells, which constants belong to which function
		std::vector<double> 			_const; 		// Values for the constants in the functions
		std::vector<std::string> 		_iso_parNames;
		std::vector<int> 			_iso_borders_par; 	// Tells, which isobar parameters belong to which isobar
		std::vector<std::string> 		_iso_constNames;
		std::vector<int> 			_iso_borders_const; 	// Tells, which isobar constants belong to which isobar
		std::vector<double> 			_iso_const; 		// Give the isobar constants

	// PHASE SPACE
		int					_globalPs; 		// Assume only ONE global PS and   //Also assume ps that does not depend on the isobar mass
		std::vector<int> 			_wavePs; 		// One PS for each wave // Or the spin

	// INTERNAL 
		bool 					_write_out;		// Flag to switch on the text_output
		std::ofstream*				_outStream;		// Stream for the text output
		bool 					_has_isobars;		// true, if de-isobarred waves are in the fit


};
#endif//WAVESET_SETS_WAVES
