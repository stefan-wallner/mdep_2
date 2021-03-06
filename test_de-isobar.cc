#include<vector>
#include<complex>

#include"anchor_t.h"


#include"matrix_utilities.h"

int main(){
	anchor_t Chi2 = anchor_t();
	Chi2.Waveset()->add_wave();
	Chi2.Waveset()->add_wave();
	Chi2.Waveset()->add_wave();
	Chi2.Waveset()->setWaveLimits(0,0.,2.);
	Chi2.Waveset()->setWaveLimits(1,0.,2.);
	Chi2.Waveset()->setWaveLimits(2,0.,2.);
	Chi2.Waveset()->setWaveName(0,"anchor");
	Chi2.Waveset()->setWaveName(1,"wave1");
	Chi2.Waveset()->setWaveName(2,"le_non_et_non");
	Chi2.Waveset()->add_func(-1);
	Chi2.Waveset()->add_func(-1);
	Chi2.Waveset()->add_func(-1);
	Chi2.Waveset()->setFunctionName(0,"anchor_func");
	Chi2.Waveset()->setFunctionName(1,"isobar_3pi");
	Chi2.Waveset()->add_iso(10);
//	Chi2.Waveset()->add_iso(10);
	Chi2.Waveset()->set_iso_const(0,1.);
	Chi2.Waveset()->setIsobarName(0,"iso1");
//	Chi2.Waveset()->setIsobarName(1,"iso2");
	Chi2.Waveset()->add_func_to_wave(0,0);
	Chi2.Waveset()->add_funcs_to_wave(1,1,0);	
//	Chi2.Waveset()->add_funcs_to_wave(1,1,1);	
	Chi2.Waveset()->add_func_to_wave(2,2);
	std::vector<double> binning;
	std::vector<double> iso_binning;
	for (int i=0;i<3;i++){
		binning.push_back(1.0*i);
	};
//	binning.push_back(3.);
	iso_binning.push_back(0.5);
	iso_binning.push_back(1.5);
	iso_binning.push_back(2.5);
	iso_binning.push_back(3.5);
	iso_binning.push_back(4.5);
	Chi2.Waveset()->setBinning(binning);
	Chi2.Waveset()->add_isobar_binning(iso_binning);
	Chi2.Waveset()->setWaveIsobarBinning(1,0);
	std::vector<double> data;
	std::vector<std::vector<double> > coma;
	for (size_t i=0;i<2*Chi2.Waveset()->nPoints()-1;i++){
		data.push_back(1.);
		std::vector<double> comaLine;
		for (size_t j=0;j<2*Chi2.Waveset()->nPoints()-1;j++){
			if(i==j){
				comaLine.push_back(1.);
			}else{
				comaLine.push_back(0.);
			};
		};
		coma.push_back(comaLine);
	};
	double data_c[] = {1.,1.,0.,4.,0.,9.,0.,16.,0.,1.,0.};
	for (int i=0;i<11;i++){
		data[i] = data_c[i];
	};
	for (size_t bin=0;bin<Chi2.Waveset()->nBins();bin++){
		if (not Chi2.set_data(0,bin,data)){
			std::cout<<"data wrong"<<std::endl;
		};
		if (not Chi2.set_coma(0,bin,coma)){
			std::cout<<"coma wrong"<<std::endl;
		};
	};
	Chi2.printStatus();

//// Start here with the evaluation
	std::vector<std::complex<double> > cpl(Chi2.Waveset()->nFtw(),std::complex<double>(1.,0.));
	std::vector<double> par;
	std::vector<std::complex<double> > bra;
	std::vector<double> isopar;
	isopar.push_back(0.);
	isopar.push_back(0.);
	isopar.push_back(1.5);
	isopar.push_back(0.);
	isopar.push_back(0.);

	cpl[1] = std::complex<double>(0.66666666666666666666,0.);
	cpl[2] = std::complex<double>(1.,0.);

	std::cout<<"EVAL"<<std::endl;
	std::cout<<"Eval: "<<Chi2.EvalCP(&cpl[0],&par[0],&isopar[0])<<std::endl;

	std::vector<std::complex<double> > startcpl(1,std::complex<double>(1.,0));

	std::cout<<"beforeGET"<<std::endl;
	AandB<double> AB = Chi2.get_AB(0,&startcpl[0],&par[0], &isopar[0]);
	std::cout<<"afterGET"<<std::endl;
	print_matrix(AB.A);
	std::vector<double> inver = cholesky::cholesky_solve(AB.A,AB.B);
	print_vector(AB.B);
	std::cout<<std::endl;
	for (size_t i=0;i<inver.size();i++){
		inver[i]*=-.5;
	};
	print_vector(inver);
	print_vector(Chi2.getMinimumCplBra(0,&bra[0],&startcpl[0],&par[0], &isopar[0]));






return 0;
};
