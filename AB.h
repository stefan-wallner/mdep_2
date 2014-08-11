#ifndef ABABABABABABABA_POP_KEK
#define ABABABABABABABA_POP_KEK
#include<Eigen/Dense>
/*
Object, that contains a n-vector and a nxn-matrix
*/
struct AandB{
		AandB(int dim);
		Eigen::MatrixXd A;
		Eigen::VectorXd B;	
};
AandB::AandB(int dim){
	A = Eigen::MatrixXd::Zero(dim,dim);
	B = Eigen::VectorXd::Zero(dim);
};
#endif//ABABABABABABABA_POP_KEK
