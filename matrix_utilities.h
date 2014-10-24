#ifndef MATRIX_UTILITIES
#define MATRIX_UTILITIES
#include<vector>
#include<iostream>

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





template<typename T>
std::vector<std::vector<T> > invert33(std::vector<std::vector<T> > mat){
	std::vector<std::vector<T> > inverse = std::vector<std::vector<T> >(3,std::vector<T>(3)); //Initialize 3x3 matrix
	T det = mat[0][0]*mat[1][1]*mat[2][2] + mat[0][1]*mat[1][2]*mat[2][0] + mat[1][0]*mat[2][1]*mat[0][2] - mat[0][2]*mat[1][1]*mat[2][0] - mat[0][1]*mat[1][0]*mat[2][2] - mat[0][0]*mat[1][2]*mat[2][1];
	if(det==0.){
		std::cerr<<"Error: invert33.h: Determinant is zero"<<std::endl;
		return inverse;
	};
	inverse[0][0] = ( mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1])/det;
	inverse[0][1] = (-mat[0][1]*mat[2][2] + mat[0][2]*mat[2][1])/det;
	inverse[0][2] = ( mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1])/det;

	inverse[1][0] = (-mat[1][0]*mat[2][2] + mat[1][2]*mat[2][0])/det;
	inverse[1][1] = ( mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0])/det;
	inverse[1][2] = (-mat[0][0]*mat[1][2] + mat[0][2]*mat[1][0])/det;

	inverse[2][0] = ( mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0])/det;
	inverse[2][1] = (-mat[0][0]*mat[2][1] + mat[0][1]*mat[2][0])/det;
	inverse[2][2] = ( mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0])/det;

	return inverse;
};

// Prints a matrix
template<typename T>
void print_matrix(const std::vector<std::vector<T> > &mat){
	std::cout<<"[";
	for (unsigned int i =0;i<mat.size();i++){
		std::cout<<"[    ";
		for (unsigned int j=0;j<mat[i].size();j++){
			std::cout << mat[i][j];
			if (j != mat[i].size()-1){
				std::cout<<",     ";
			};
		};
		if (i != mat.size()-1){
			std::cout<<"],"<<std::endl;
		}else{
			std::cout<<"]]"<<std::endl;
		};
	};
};
// Prints a vector
template<typename T>
void print_vector(const std::vector<T> &in){
	std::cout << "[";
	if (in.size()>0){
		std::cout << in[0];
		for (unsigned int i=1;i<in.size();i++){
			std::cout << ", " << in[i];
		};
	};
	std::cout<<"]"<<std::endl;
};

// Gets the coefficients of a parabola from three points (x[i],y[i]) i in [0,1,2])
template<typename T> 
std::vector<T> get_abc(std::vector<T> &x, std::vector<T> &y){
	std::vector<std::vector<T> > matrix;
	std::vector<T> squares;
	std::vector<T> linears;
	std::vector<T> ones;
	for (int i=0; i<3;i++){
		squares.push_back(x[i]*x[i]);
		linears.push_back(x[i]);
		ones.push_back(1.);
	};
	matrix.push_back(squares);
	matrix.push_back(linears);
	matrix.push_back(ones);
	std::vector<std::vector<T> > inverse = invert33(matrix);
	std::vector<T> abc;
	for (int i=0;i<3;i++){
		T val = 0.;
		for (int j=0;j<3;j++){
			val+=inverse[j][i]*y[j];
		};
		abc.push_back(val);
	};
	return abc;
};

// At the moment no error handling is included, 
// The sizes of the matrices, as well as their symmetricity and their positive definiteness have to be ensured
// Also seems to work for complex matrices 

namespace cholesky{
	template<typename xdouble>
	void print_vector(std::vector<xdouble> &x){
		int dim = x.size();
		for (int i=0;i<dim;i++){
			std::cout<<x[i]<<"\t";
		};
		std::cout<<std::endl<<std::endl;
	};

	template<typename xdouble>
	void print_matrix(std::vector<std::vector<xdouble> > &A){
		int dim = A.size();
		for (int i =0;i<dim;i++){
			for (int j=0;j<dim;j++){
				std::cout<<A[i][j]<<"\t";
			};
			std::cout<<std::endl;
		};
		std::cout<<std::endl;
	};

	template<typename xdouble> // C = A*B // c_{ij} = sum_k a_{ik}b_{kj}
	std::vector<std::vector<xdouble> > dot(std::vector<std::vector<xdouble> > &A, std::vector<std::vector<xdouble> > &B){
		int dim = A.size();
		std::vector<std::vector<xdouble> > C(dim,std::vector<xdouble>(dim,0.));
		for (int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){
				for (int k=0;k<dim;k++){
					C[i][k] += A[i][j]*B[j][k];
				};
			};
		};
		return C;
	};

	template<typename xdouble> // y = A*x // y_i = sum_j A_{ij}x_j
	std::vector<xdouble> dot(std::vector<std::vector<xdouble> > &A, std::vector<xdouble> &x){
		int dim = x.size();
		std::vector<xdouble> y = std::vector<xdouble>(dim,0.);
		for (int i =0;i<dim;i++){
			for (int j=0;j<dim;j++){
				y[i]+=A[i][j]*x[j];
			};
		};
		return y;
	};

	template<typename xdouble> // T = A^T
	std::vector<std::vector<xdouble> > transpose(std::vector<std::vector<xdouble> > &A){
		int dim = A.size();
		std::vector<std::vector<xdouble> > T(dim,std::vector<xdouble>(dim));
		for (int i=0;i<dim;i++){
			for (int j=0;j<dim;j++){
				T[i][j] = A[j][i];
			};
		};
		return T;
	};

	template<typename xdouble> // A = G*G^T // This destroys A!!!
	std::vector<std::vector<xdouble> > cholesky_decompose(std::vector<std::vector<xdouble> > &A){
		int dim = A.size();
		std::vector<std::vector<xdouble> > G(dim,std::vector<xdouble>(dim,0.));
		for (int j=0;j<dim;j++){
			G[j][j] = sqrt(A[j][j]);
			for (int i=j+1;i<dim;i++){
				G[i][j] = A[i][j]/G[j][j];
				for (int k=j+1;k<i+1;k++){
					A[i][k]-=G[i][j]*G[k][j];
				};
			};
		};
		return G;
	};

	template<typename xdouble> // Solves A*x = b, assuming A is an upper triangular matrix
	std::vector<xdouble> upper_diag_solve(std::vector<std::vector<xdouble> > &A, std::vector<xdouble> &b){
		int dim = b.size();
		std::vector<xdouble> x(dim);
		for (int i=dim-1;i >=0;i--){
			x[i] = b[i];
			for (int j=i+1;j<dim;j++){
				x[i]-=A[i][j]*x[j];
			};
			x[i]/=A[i][i];
		};
		return x;
	};

	template<typename xdouble> // Solves A*x = b, assuming A is an lower triangular matrix
	std::vector<xdouble> lower_diag_solve(std::vector<std::vector<xdouble> > &A, std::vector<xdouble> &b){
		int dim = b.size();
		std::vector<xdouble> x(dim);
		for (int i=0;i<dim;i++){
			x[i] = b[i];
			for (int j=0;j<i;j++){
				x[i]-=A[i][j]*x[j];
			};
			x[i]/=A[i][i];
		};
		return x;
	};

	template<typename xdouble> // Solves A*x=B assuming A is symmetrix and positive definite // This destroys A!!!
	std::vector<xdouble> cholesky_solve(std::vector<std::vector<xdouble> > &A, std::vector<xdouble> &b){
		std::vector<std::vector<xdouble> > G = cholesky_decompose(A);
		std::vector<std::vector<xdouble> > GT= transpose(G);
		std::vector<xdouble> y = lower_diag_solve(G,b);
		std::vector<xdouble> x = upper_diag_solve(GT,y);
		return x;
	};

	template<typename xdouble> // Checks if A is symmetric
	bool is_symmetric(std::vector<std::vector<xdouble> > &A){
		int dim = A.size();
		for (int i=0;i<dim;i++){
			for (int j=0;i<i;j++){
				if(A[i][j] != A[j][i]){
					return false;
				};
			};
		};
		return true;
	};

};
#endif//MATRIX_UTILITIES
