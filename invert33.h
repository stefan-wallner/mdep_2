#ifndef INVERT33
#define INVERT33
#include<vector>
#include<iostream>
// Calculates the invese of a 3x3 matrix
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

#endif//INVERT33
