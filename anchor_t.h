#ifndef ANCHOR_ROHCNA
#define ANCHOR_ROHCNA
#include"waveset.h"
#include<complex>
#include<vector>
#include<string>
#include"AB.h"


class anchor_t : public waveset {
	public:
		std::vector<std::complex<double> > get_branchings(std::vector<std::complex<double> > &cpl,std::vector<double> &par, std::vector<double> &iso_par);

		std::vector<std::vector<double> > getPlots(int tbin, std::vector<std::complex<double> > &cpl, std::vector<double> &par);

		anchor_t();

		template<typename xdouble>
		xdouble EvalCP(std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par );

		template<typename xdouble>
		xdouble EvalBranch(std::vector<std::complex<xdouble> >&branch, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par);

		template<typename xdouble>
		xdouble EvalAutoCpl(std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par, std::vector<xdouble> &iso_par);

		template<typename xdouble>
		xdouble EvalAutoCplBranch(std::vector<std::complex<xdouble> >&bra, std::vector<std::complex<xdouble> >&cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par);

		template<typename xdouble> 
		xdouble EvalAutoCplTbin(int tbin, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par);

		template<typename xdouble>
		xdouble EvalTbin(int tbin,std::vector<std::complex<xdouble> > &cpl,std::vector<xdouble> &par, std::vector<xdouble> &iso_par);

		template<typename xdouble>
		xdouble EvalBin(int tbin, int bin, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par, std::vector<std::vector<std::complex<xdouble> > > &iso_eval);

		template<typename xdouble>
		std::vector<xdouble> delta(int tbin, int bin,double mass, std::vector<std::complex<xdouble> > &cpl, std::vector<xdouble> &par, std::vector<std::vector<std::complex<xdouble> > > &iso_eval);

		std::vector<double> get_data(int tbin, int bin);
		std::vector<std::vector<double> > get_coma(int tbin, int bin);


		void loadData(int tbin, const char* dataFile);
		void loadComa(int tbin, const char* comaFile);
		void nullify();
		void conjugate();

		int getNtotAnc();
		int getNanc();

		std::string className();

		void update_n_cpls();
		template<typename xdouble>
		AandB<xdouble> get_AB(int tbin,std::vector<std::complex<xdouble> > &anchor_cpl, std::vector<xdouble> &par);

		template<typename xdouble>
		AandB<xdouble> get_AB_iso(int tbin, std::vector<std::complex<xdouble> > &anchor_cpl,std::vector<xdouble> &par, std::vector<xdouble> &iso_par); // For de-isobarred

		template<typename xdouble>
		std::vector<std::complex<xdouble> > getMinimumCpl(int tbin,std::vector<std::complex<xdouble> > &anchor_cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par);

		template<typename xdouble>
		std::vector<std::complex<xdouble> > getMinimumCplBra(int tbin, std::vector<std::complex<xdouble> > &branch, std::vector<std::complex<xdouble> > &anchor_cpl, std::vector<xdouble> &par, std::vector<xdouble> &iso_par);

		bool set_data(int tbin, int bin, std::vector<double> data);
		bool set_coma(int tbin, int bin, std::vector<std::vector<double> > coma);
		void printStatus();
		void set_is_ampl(bool is_ampl);
		void setTbinning(std::vector<double> binning);

	protected:
		int _nBrCplAnc;// Number of couplings with branchings in the anchor wave
		bool _is_ampl; // Enable Amplitunde fitting, needs to be tested.
		std::vector<std::vector<std::vector<double> > > _data; // data 
		std::vector<std::vector<std::vector<std::vector<double> > > > _coma; // covariance matrix
};


#endif//ANCHOR_ROHCNA
