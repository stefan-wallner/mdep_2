#ifndef ABABABABABABABA_POP_KEK
#define ABABABABABABABA_POP_KEK
#include<vector>
/*
Object, that contains a n-vector and a nxn-matrix
*/
template<typename xdouble>
struct AandB{
		AandB(int dim);
		std::vector<std::vector<xdouble> > A;
		std::vector<xdouble> B;	
};

template<typename xdouble>
AandB<xdouble>::AandB(int dim){
	A = std::vector<std::vector<xdouble> >(dim,std::vector<xdouble>(dim,0.));
	B = std::vector<xdouble>(dim,0.);
};
#endif//ABABABABABABABA_POP_KEK
