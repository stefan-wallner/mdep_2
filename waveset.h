#ifndef WAVESET_SETS_WAVES
#define WAVESET_SETS_WAVES
#include<vector>
#include<complex>
#include<fstream>
#include<iostream>

#ifdef USE_YAML
#include "yaml-cpp/yaml.h"
#endif//USE_YAML

#include"amplitude.h"

class waveset {
	public:

	// CONSTRUCTOR
		waveset();
#ifdef USE_YAML
		waveset(
								std::string 						card);
#endif//USE_YAML
	// FUNCTIONS AND AMPLITUDES

		template<typename xdouble>
		std::vector<std::complex<xdouble> > amps(
								double 							m,
								const std::complex<xdouble> 				*cpl,
								const xdouble	 					*par,
								std::vector<std::vector<std::complex<xdouble> > > 	&funcEvals2pi) 	const;

		template<typename xdouble>
		std::vector<std::complex<xdouble> > funcs(
								double 							m,
								const xdouble 						*par) 		const;

		template<typename xdouble>
		std::vector<std::vector<std::complex<xdouble> > > iso_funcs(
								const xdouble 						*par) 		const;

		std::vector<double> phase_space(
								double 							m) 		const;

	// SET UP WAVESET
		// // WAVES AND FUNCTIONS
		size_t 				add_wave();
		size_t 				add_func(int i);
		size_t 				add_iso(int i);

		// // SET AMPLITUDE DEFINITIONS
		void 				add_func_to_wave(int wave, int func);
		void 				add_funcs_to_wave(int wave, int func, int func2);
		void 				couple_funcs(int i, int j);

		// // SETTINGS FOR DIFFERENT WAVES
		void 				setWaveLimits(int i, double lower, double upper);
		void 				setWaveSpin(int i, int L);
		void 				setWaveIsobarSpin(int wave, int L);
		void 				setGlobalPhaseSpace(int i);
		void 				setWavePhaseSpace(int i, int ps);
		void 				setWaveIsobarBinning(int wave, int binning);
		void 				setConst(int i,double con);
		void 				set_iso_const(int con, double value);
		void 				add_isobar_binning(std::vector<double> binning);

		// // NAMES
		void 				setWaveName(int i, std::string name);
		void 				setFunctionName(int i, std::string name);
		void 				setParameterName(int i, std::string name);
		void 				setConstantName(int i, std::string name);
		void 				setIsobarName(int i,std::string name);
		void 				setIsoParName(int i,std::string name);
		void 				setIsoConstName(int i,std::string name);

		// // BINNING
		void 				setBinning(std::vector<double> binning);
		void 				setTbinning(std::vector<double> binning);
		void 				setEvalTbin(int i, bool flag);
		void 				setMinBin(int in);
		void 				setMaxBin(int in);


#ifdef USE_YAML
		// // YAML LOADER
		bool				loadGlobalPhaseSpace(YAML::Node &waveset);
		std::map<std::string,int>	loadFunctions(YAML::Node &waveset, YAML::Node &param);
		bool 				loadWaves(YAML::Node &waveset, YAML::Node &defs);
		void				loadFtw(YAML::Node &waveset, std::map<std::string,int> &fMap);
		void				loadBranchings(YAML::Node &waveset);
		void				loadBinnings(YAML::Node &waveset);
		void				loadIsoBinnings(YAML::Node &waveset, YAML::Node &iso_binnings);
#endif//USE_YAML

	// GETTERS

		// // ONE LINE GETTERS
		bool 				write_out		()const		{return _write_out;};			
		bool				has_isobars		()const		{return _has_isobars;};				
		int				minBin			()const		{return _minBin;};			
		int				maxBin			()const		{return _maxBin;};			
		int				nBrCpl			()const		{return _nBrCpl;};			
		// // SIMPLE OVER ALL NUMBERS
		int 				nPoints			()const		{return _nPoints;};
		int 				nFtw			()const		{return _nFtw;};
		int 				nTbin			()const		{return _nTbin;};
		int 				nBins			()const		{return _nBins;};
		int 				getNtot			()const		{return 2*getNcpl() + getNpar() + 2*nBranch() + getNiso();};
		int 				getNcpl			()const		{return _nBrCpl * _nTbin;};
		int 				nBranch			()const		{return _nBranch;};
		double				get_m			(int bin)const	{return (_binning[bin]+_binning[bin+1])/2.;};
		std::ofstream*			outStream		()const		{return _outStream;};
		const std::vector<bool>*	eval_tbin		()const		{return &_eval_tbin;};

		const std::vector<int>* 	n_branch		()const		{return &_n_branch;};
		const std::vector<int>*		n_cpls			()const		{return &_n_cpls;};
		const std::vector<int>*		borders_waves		()const		{return &_borders_waves;};
		const std::vector<int>*		point_to_wave		()const		{return &_point_to_wave;};
		const std::vector<int>*		funcs_to_waves		()const		{return &_funcs_to_waves;};
		const std::vector<int>*		iso_to_waves		()const		{return &_iso_to_waves;};
		const std::vector<int>*		wave_binning_pts	()const		{return &_wave_binning_pts;};
		const std::vector<int>*		point_borders_wave	()const		{return &_point_borders_wave;};
		const std::vector<int>*		coupled			()const		{return &_coupled;};

		const std::vector<double>*	upperLims		()const		{return &_upperLims;};
		const std::vector<double>*	lowerLims		()const		{return &_lowerLims;};
		const std::vector<double>*	binning			()const		{return &_binning;};
		const std::vector<double>*	t_binning		()const		{return &_t_binning;};

		// // PROPERTIES OF THE WAVES
		int 				getNpar()				const;
		int 				getNiso()				const;
		std::string 			getWaveName(int i)			const;
		std::vector<int>	 	get_wave_functions(int wave)		const;
		std::vector<int>	 	get_wave_pars(int wave)			const;
		std::vector<int>	 	get_wave_const(int wave)		const;
		std::vector<int>	 	get_wave_isobars(int wave)		const;
		std::vector<int>	 	get_wave_iso_pars(int wave)		const;
		std::vector<int>	 	get_wave_iso_const(int wave)		const;

		// // PROPERTIES OF THE FUNCTIONS
		std::string 			getFunctionName(int i)			const;
		std::vector<int>	 	get_nPars()				const;
		std::vector<int>	 	get_nConst()				const;
		std::vector<int>	 	get_function_pars(int func)		const;
		std::vector<int>	 	get_function_const(int func)		const;
		std::vector<int>	 	get_function_waves(int func)		const;

		// // PROPERTIES OF THE ISOBARS
		std::string 			getIsobarName(int i)			const;
		std::vector<int>	 	get_nParsIso()				const;
		std::vector<int>	 	get_nConstIso()				const;
		std::vector<int>	 	get_isobar_pars(int func)		const;
		std::vector<int>	 	get_isobar_const(int func)		const;
		std::vector<int>	 	get_isobar_waves(int func)		const;

		// // PROPERTIES OF THE PARAMETERS
		std::string 			getParameterName(int i)			const;
		std::string 			getConstantName(int i)			const;
		std::string 			getIsoParName(int i)			const;
		std::string 			getIsoConstName(int i)			const;

		// // POPERTIES OF THE BINNING & BRANCHING
		int 				get_bin(double mass)			const;
		std::vector<int>	 	getFirstBranch()			const;

	// UPDATERS
		void 				updateNftw();
		void 				updateNpoints();
		void 				updateFuncLims();
		void 				updateFuncSpin();
		void 				updateIsobar();
		void 				update_min_max_bin();
		void 				handle_branchings(int wave, int func);
		void 				update_n_cpls();
		void 				update_n_branch();
		void 				updateTprime(int tbin);

	// INFO FUNCTIONS
		std::string 			className		()const		{return "waveset";};
		bool 				checkConsistency()			const;
		void 				printStatus()				const;
		void 				printParameters()			const;
		void 				open_output(std::string filename ="chi2log.dat");
		void 				close_output();

	protected:
		double					_t_prime;		// Four momentum transfer
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
		std::vector<amplitude*>			_amp_funcs;		// Defined amplitude functions
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
		std::vector<std::string> 		_iso_parNames;
		std::vector<int> 			_iso_borders_par; 	// Tells, which isobar parameters belong to which isobar
		std::vector<std::string> 		_iso_constNames;
		std::vector<int> 			_iso_borders_const; 	// Tells, which isobar constants belong to which isobar
		std::vector<double> 			_iso_const; 		// Give the isobar constants

	// PHASE SPACE
		int					_globalPs; 		// Assume only ONE global PS and   //Also assume ps that does not depend on the isobar mass
		std::vector<int> 			_wavePs; 		// One PS for each wave // Or the spin

	// BINNING
		int 					_nBins; 		// Number of bins
		std::vector<double> 			_binning;   		// Binning in m3pi
		int 					_minBin; 		// Minimum bin used by any wave
		int 					_maxBin; 		// Maximum bin used by any wave
		double 					_mMin; 			// minimum mass (m3pi)
		double 					_mMax; 			// maximum mass (m3pi)
		int 					_nTbin;  		// number of t' bins
		std::vector<double> 			_t_binning; 		// Binning in t'
		std::vector<bool> 			_eval_tbin; 		// switch on/off single t' bins

	// BRANCHING
		int 					_nBranch;  		// Number of branchings
		int 					_nBrCpl;   		// Number of couplings with branchings
		std::vector<int> 			_coupled;  		// Encodes coupled functions
		std::vector<int>			_n_branch; 		// Number of branching for wave/function
		std::vector<int> 			_n_cpls;   		// Number of coupling for wave/function

	// INTERNAL
		bool 					_write_out;		// Flag to switch on the text_output
		std::ofstream*				_outStream;		// Stream for the text output
		bool 					_has_isobars;		// true, if de-isobarred waves are in the fit
};
#endif//WAVESET_SETS_WAVES
