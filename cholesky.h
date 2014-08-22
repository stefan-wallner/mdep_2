#ifndef CHOLESKY_SUPER_DUPER
#define CHOLESKY_SUPER_DUPER
#include<vector>
#include<math.h>

// At the moment no error handling is included, 
// The sizes of the matrices, as well as their symmetricity and their positive definiteness have to be ensured

namespace cholesky{
	template<typename xdouble>
	void print_vector(std::vector<xdouble> x){
		int dim = x.size();
		for (int i=0;i<dim;i++){
			std::cout<<x[i]<<"\t";
		};
		std::cout<<std::endl<<std::endl;
	};

	template<typename xdouble>
	void print_matrix(std::vector<std::vector<xdouble> > A){
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
	std::vector<std::vector<xdouble> > dot(std::vector<std::vector<xdouble> > A, std::vector<std::vector<xdouble> > B){
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
	std::vector<xdouble> dot(std::vector<std::vector<xdouble> > A, std::vector<xdouble> x){
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
	std::vector<std::vector<xdouble> > transpose(std::vector<std::vector<xdouble> > A){
		int dim = A.size();
		std::vector<std::vector<xdouble> > T(dim,std::vector<xdouble>(dim));
		for (int i=0;i<dim;i++){
			for (int j=0;j<dim;j++){
				T[i][j] = A[j][i];
			};
		};
		return T;
	};

	template<typename xdouble> // A = G*G^T
	std::vector<std::vector<xdouble> > cholesky_decompose(std::vector<std::vector<xdouble> > A){
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
	std::vector<xdouble> upper_diag_solve(std::vector<std::vector<xdouble> > A, std::vector<xdouble> b){
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
	std::vector<xdouble> lower_diag_solve(std::vector<std::vector<xdouble> > A, std::vector<xdouble> b){
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

	template<typename xdouble> // Solves A*x=B assuming A is symmetrix and positive definite
	std::vector<xdouble> cholesky_solve(std::vector<std::vector<xdouble> > A, std::vector<xdouble> b){
		if (not is_symmetric(A)){
			std::cout<<"Try to fool me, bro?"<<std::endl;
		};
		std::vector<std::vector<xdouble> > G = cholesky_decompose(A);
		std::vector<std::vector<xdouble> > GT= transpose(G);
		std::vector<xdouble> y = lower_diag_solve(G,b);
		std::vector<xdouble> x = upper_diag_solve(GT,y);
		return x;
	};

	template<typename xdouble> // Checks if A is symmetric
	bool is_symmetric(std::vector<std::vector<xdouble> > A){
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


#endif//CHOLESKY_SUPER_DUPER
