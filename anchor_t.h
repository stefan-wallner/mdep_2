#ifndef ANCHOR_ROHCNA
#define ANCHOR_ROHCNA
#include"chi2_2d.h"
#include<complex>
#include<vector>
#include<string>
#include"AB.h"


class anchor_t : public chi2_2d {
	public:
		anchor_t();
		
		template<typename xdouble>
		xdouble EvalCP(std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par);

		template<typename xdouble>
		xdouble EvalBranch(std::vector<std::complex<xdouble> >&branch, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par);

		double EvalAutoCpl(std::vector<std::complex<double> > &cpl,std::vector<double> &par);
		double EvalAutoCplBranch(std::vector<std::complex<double> >&bra, std::vector<std::complex<double> >&cpl, std::vector<double> &par);
		double EvalAutoCplTbin(int tbin, std::vector<std::complex<double> > &cpl, std::vector<double> &par);

		template<typename xdouble>
		xdouble EvalTbin(int tbin,std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par);

		template<typename xdouble>
		xdouble EvalBin(int tbin, int bin, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par);


		template<typename xdouble>
		std::vector<xdouble> delta(int tbin, int bin,double mass, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par);

		int getNtot();
		int getNtotAnc();
		int getNcpl();
		int getNanc();
		int getNpar();
		int getNbra();
		int getNiso();

		int getNtBin();

		void setBinning(std::vector<double> binning);
		void setTbinning(std::vector<double> binning);
		void loadData(int tbin, const char* dataFile);
		void loadComa(int tbin, const char* comaFile);
		void nullify();
		void conjugate();
		
		void add_func(int i, bool ist_t_dep = false);
		void setConst(int i,double con);
		void add_func_to_wave(int wave, int func);
		void couple_funcs(int i, int j);
		void setWaveLimits(int i, double lower, double upper);

		std::string className();

		std::vector<std::vector<double> > getPlots(int tbin, std::vector<std::complex<double> > &cpl, std::vector<double> &par);


		AandB get_AB(int tbin,std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par);
		std::vector<std::complex<double> > getMinimumCpl(int tbin,std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par);
		std::vector<std::complex<double> > getMinimumCplBra(int tbin, std::vector<std::complex<double> > &branch, std::vector<std::complex<double> > &anchor_cpl, std::vector<double> &par);
		void updateTprime(int tbin);

		int get_bin(double mass);
		std::vector<double> get_data(int tbin, int bin);
		std::vector<std::vector<double> > get_coma(int tbin, int bin);

		std::vector<std::complex<double> > get_branchings(std::vector<std::complex<double> > cpl,std::vector<double> par);

		void update_n_cpls();
		void update_n_branch();
		void update_min_max_bin();

		void printStatus();
		std::vector<int> getFirstBranch();

		void setEvalTbin(int i, bool flag);
		void set_is_ampl(bool is_ampl);
	protected:
		bool _is_ampl; // Enable Amplitunde fitting, needs to be tested.
		int _nBins; // Number of bins
		int _minBin; // Minimum bin used by any wave
		int _maxBin; // Maximum bin used by any wave
		int _nTbin;  // number of t' bins
		double _mMin; // minimum mass (m3pi)
		double _mMax; // maximum mass (m3pi)
		std::vector<int> _const_is_t; // List of constants, that are t' actually (will then be set automatically)
		std::vector<double> _t_binning; // Binning in t'
		std::vector<double> _binning;   // Binning in m3pi
		std::vector<std::vector<std::vector<double> > > _data; // data 
		std::vector<std::vector<std::vector<std::vector<double> > > > _coma; // covariance matrix

		int _nBranch;  // Number of branchings
		int _nBrCpl;   // Number of couplings with branchings
		int _nBrCplAnc;// Number of couplings with branchings in the anchor wave
		std::vector<int> _coupled;  // Encodes coupled functions
		std::vector<int> _n_branch; // Number of branching for wave/function
		std::vector<int> _n_cpls;   // Number of coupling for wave/function
					    // Map of how to treat the nonAnchor couplings in the analytic calculation.
		std::vector<bool> _eval_tbin; // switch on/off single t' bins


};


#endif//ANCHOR_ROHCNA
