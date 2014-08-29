#ifndef CHI2_LOLOLO2D
#define CHI2_LOLOLO2D
#include"chi2.h"
#include<vector>
#include<complex>

class chi2_2d : public chi2{
	public:
		// CONSTRUCTORS
		chi2_2d();
		// EVALUATE THE FUNCTION

		template<typename xdouble>
		std::vector<std::complex<xdouble> > amps(double m,std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par, std::vector<std::vector<std::complex<xdouble> > > &funcEvals2pi);

		template<typename xdouble>
		std::vector<std::vector<std::complex<xdouble> > > iso_funcs(std::vector<xdouble> &par);

//		std::vector<double> phase_space(double m);
		// MANAGE FUNCTIONS
		std::string className();
		void add_wave();
		void add_func_to_wave(int wave, int func);
		void add_funcs_to_wave(int wave, int func, int func2);


		// HELPERS
		bool checkConsistency();
		void printStatus();
		void printParameters();

		// METHODS FOR INTERFACING WITH THE 'OLD' PARAMETERS [RE,IM,PAR,...,RE,IM,PAR...]
////////////// 2d specials ////////////////////
		int getNpoints();

		void updateIsobar();

		void add_isobar_binning(std::vector<double> binning);
		void setWaveIsobarSpin(int wave, int L);
		void setWaveIsobarBinning(int wave, int binning);
		void add_iso(int i);
		void updateNpoints();

		std::vector<int> get_wave_isobars(int wave);
		std::vector<int> get_wave_iso_pars(int wave);
		std::vector<int> get_wave_iso_const(int wave);
		std::vector<int> get_isobar_pars(int func);
		std::vector<int> get_isobar_const(int func);
		std::vector<int> get_isobar_waves(int func);

		std::vector<int> get_nParsIso();
		std::vector<int> get_nConstIso();

		void set_iso_const(int con, double value);

		void setIsobarName(int i,std::string name);
		void setIsoParName(int i,std::string name);
		void setIsoConstName(int i,std::string name);
		std::string getIsobarName(int i);
		std::string getIsoParName(int i);
		std::string getIsoConstName(int i);

	protected:

		int _nPoints; 					// Number of amplitudes in the end (+1 for each wave, +1 for each isobar step)
		bool _has_isobars;				// true, if de-isobarred waves are in the fit

		int _nIso; 					// Number of isobars
		int _maxNparsIso; 				// Maximum number of parameters for one isobar
		int _maxBinIso; 				// Maximum number of bins for one isobar
		std::vector<int> _isos; 			// Isobar paramterizations
		std::vector<int> _iso_to_waves; 		// Couples isobars to waves (-1 -> no isobar paramtrization -> Standard mdep fit)
		std::vector<int> _iso_borders_par; 		// Tells, which isobar parameters belong to which isobar
		std::vector<int> _iso_borders_const; 		// Tells, which isobar constants belong to which isobar
		std::vector<int> _L_iso; 			// Gives the spin of the isobar for the single waves
		std::vector<int> _L_iso_func; 			// Gives the spin of the the single isobars
		std::vector<double> _iso_const; 		// Give the isobar constants
		std::vector<int> _iso_binning_pts; 		// Number of isobar bins for each wave (no isobar -> -1)
		std::vector<int> _wave_binning_pts; 		// Number of bins for the isobars
		std::vector<int> _point_to_wave;		// Tells, which points belong to which wave
		std::vector<int> _iso_n_binning; 		// Tells, which isobar binning to use for the waves
		std::vector<int> _wave_n_binning; 		// Tells, which isobar binning for the single isobars
		std::vector<std::vector<double> > _iso_binnings;// Gives different isobar binnings
		std::vector<int> _point_borders_wave;		// Tells, where the isobar points for a single wave start


		// Names
		std::vector<std::string> _iso_funcNames;
		std::vector<std::string> _iso_parNames;
		std::vector<std::string> _iso_constNames;
};
#endif//CHI2_LOLOLO2D
