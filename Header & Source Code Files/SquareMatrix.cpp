/*
 * SquareMatrix.cpp
 *
 *  Created on: 21.12.17
 *      Author: Edward James - Medical Imaging CDT MRes, 2017-2018. 
 */

// Source file for class SquareMatrix, to define implementations for functions declared in SquareMatrix.h

#include "SquareMatrix.h"
#include <cstdlib>
#include <cmath>  //for abs function

using namespace std;

// Part 2. 

// 1. ConstructorsTesting

//default constructor	
SquareMatrix::SquareMatrix () : Matrix () {
}

//standard custom constructor   
SquareMatrix::SquareMatrix (const int& dim) : Matrix(dim,dim) {
}

//default copy constructor #1	
SquareMatrix::SquareMatrix (const SquareMatrix& input) : Matrix(input) {   //returns SquareMatrix copy of SquareMatrix
}

//default copy constructor #2	
SquareMatrix::SquareMatrix (const Matrix& input) : Matrix(input) {   //can return SquareMatrix copy of Matrix if Matrix is square
	//error check to see if attempting to return square copy of non-square matrix ...
	checkSquare();	
}

void SquareMatrix::checkSquare () const {
	if (this->noOfRows != this->noOfCols) { 
		cerr << endl;		
		cerr << "Error:  In default copy constructor #2 ..." << endl;
		cerr << "	Cannot create a square matrix with unequal dimensions." << endl;
		cerr << "        Exiting program ... " << endl;
		cerr << endl;
		exit(1);		
	}
}

//static zeros constructor 	
SquareMatrix SquareMatrix::Zeros(const int& dim) {
	
	//return an instance of type SquareMatrix invoked by static zeros function in class Matrix of appropriate dimensions
	return Matrix::Zeros(dim,dim);
}

//static ones constructor	
SquareMatrix SquareMatrix::Ones(const int& dim) {
	
	//return an instance of type SquareMatrix invoked by static ones function in class Matrix of appropriate dimensions
	return Matrix::Ones(dim,dim);
}

//static identity matrix constructor 
SquareMatrix SquareMatrix::Eye(const int& dim) {
	
	// ceate a SquareMatrix object of zeros using the standard custom constructor
	SquareMatrix Eye(dim);	
	
	//loop through the row indices of Eye
	for (int i = 0; i < dim; ++i) {
		//looping through the column indices of Eye
		for (int j = 0; j < dim; ++j) {
			int index = i + j*dim; //to index Eye data   (cannot use GetIndex(i,j) here)
			if (i==j) {
				Eye.data[index] = 1; //assign '1' to entries on main diagonal
			}
			else {		
				Eye.data[index] = 0; //assign '0' to other entries
			}			
		}
	}	
	
	// return the object
	return Eye;
}

//for constructing specified test square matrices for general testing
SquareMatrix SquareMatrix::Test (const double* data_rhs, const int& dim) { 
	
	return Matrix::Test (data_rhs, dim, dim);
}

//destructor 
SquareMatrix::~SquareMatrix () {
	//will always invoke base class destructor, therefore do not delete dynamic memory here
}

// 2. ToeplitzTesting & ToeplitzTestingHelper

//creates a new Toeplitz square matrix when only first row is specified
SquareMatrix SquareMatrix::Toeplitz (const double* const row, const int& dim) {  
	
	return Matrix::Toeplitz (row, dim, row, dim);
}

//creates a new Toeplitz square matrix when both first column and first row are specified
SquareMatrix SquareMatrix::Toeplitz (const double* const column, const double* const row, const int& dim) {

	return Matrix::Toeplitz (column, dim, row, dim);
}

// 3. TransposeTesting

//no new function declarations required.

// 4. TriangularExtractionTesting
	
//extract upper triangular part of square matrix, converts all values below the main diagonal to zero, non-static member version
SquareMatrix SquareMatrix::TriUpper () const {	

	//initialise result with invoking matrix using default copy constructor
	SquareMatrix result(*this);	

	int dim = this->noOfRows;

	// looping through the row indices of square matrix	
	for (int i = 0; i < dim; ++i) {
		//looping through the column indices of square matrix
		for (int j = 0; j < dim; ++j) {
			int index= GetIndex (i, j);  //to index entries
			if (i>j) {			
				result.data[index] = 0;    //set values below the main diagonal to zero
			}
		}
	}

	return result;
}
	
//extract lower triangular part of square matirx, converts all values above the main diagonal to zero, non-static member version
SquareMatrix SquareMatrix::TriLower () const {	
	
	//initialise result with invoking matrix using default copy constructor
	SquareMatrix result(*this);
	
	int dim = this->noOfRows;

	// looping through the row indices of square matrix	
	for (int i = 0; i < dim; ++i) {
		//looping through the column indices of square matrix
		for (int j = 0; j < dim; ++j) {
			int index= GetIndex (i, j);	//to index entries
			if (i<j) {			
				result.data[index] = 0;	//set values above the main diagonal to zero
			}
		}
	}

	return result;
}

// 5. LUDecompositionTesting - Algorithm 2.1 - Gaussian Elimination without Pivoting

//Non-Static version read only function
SquareMatrix SquareMatrix::LUDecompositionOne (const char& output) const {	

	//initialise U with invoking matrix using default copy constructor
	SquareMatrix U(*this);	
	
	//initialise L with an identity matrix of appropriate dimension
	int dim = this->noOfRows;
	SquareMatrix L = SquareMatrix::Eye(dim);
	
	//looping through the column indices of L and U
	for (int k = 0; k < dim-1; ++k) {  
		//looping through the row indices of L and U
		for (int j = k+1; j < dim; ++j) {  	
			
			int index_jk = GetIndex(j,k);	//to index jk entries  
			int index_kk = GetIndex(k,k);   //to index kk entries
			
			//to check if diagonal entries in U (aka pivots) are zero
			if (U.data[index_kk] == 0) {
				cerr << endl;
				cerr << "Error:  In Algorithm 2.1 - Gaussian Elimination without Pivoting ..." << endl;
				cerr << "\tA zero pivot value has occurred." << endl;				
				cerr << "\tAlgorithm failed." << endl;
				cerr << endl;
				SquareMatrix error;   //use the default constructor to return the minimal version of SquareMatrix
				return error;
			}

			//if not then proceed with algorithm  (NB - no tmp data holders required here)
			else {
				//to alter L_jk value					
				L.data[index_jk] = U.data[index_jk]/U.data[index_kk];
			
				//to alter U values in row j from columns k to dim
				for (int i = k; i < dim; ++i) {		
					int index_ji = GetIndex(j,i);	//to index j,k:m entries
					int index_ki = GetIndex(k,i);	//to index k,k:m entries			
					U.data[index_ji] -=  (L.data[index_jk]*U.data[index_ki]);
				}			
			}
		}
	}

	if (output == 'L') {
		return L;
	}
	if (output == 'U') {
		return U;
	}

}

// 6. LUDecompositionTesting - Algorithm 2.2 - Gaussian Elimination with Partial Pivoting

//Non-static read only function
SquareMatrix SquareMatrix::LUDecompositionTwo (const char& output) const {	
	
	//initialise U with invoking matrix using default copy constructor
	SquareMatrix U(*this);	
	
	//initialise L and P with an identity matrix of appropriate dimension
	int dim = this->noOfRows;
	SquareMatrix L = SquareMatrix::Eye(dim);
	SquareMatrix P = SquareMatrix::Eye(dim);	

	//looping through the column indices of L, U and P
	for (int k = 0; k < dim-1; ++k) {
		
		// select i>=k to maximise magnitude of U_ik  (i.e. for column k, determine entry, on or below the diagonal of U, with the largest absolute value)
		int i = 0;
		double max = 0;		
		for (int b = k; b < dim; ++b) {			//use b instead of i here to avoid confusion
			int index_bk = GetIndex(b,k);		//to index bk entries  
			double tmp = abs(U.data[index_bk]);
			if (tmp > max) {
				max = tmp;
				i = b;
			}			
		}

		// U_ik is now the new pivot in column k
		// Swap the row (k) containing the diagonal of the column (the default pivot) with the row (i) containing the new pivot, in all 3 matrices:
		
		//interchange U_(k,k:m) and U_(i,k:m) (i.e. swap row k with row i in U matrix, but only between columns k to m)
		U.ExchangeRows(k, i, k, dim-1);

		//interchange L_(k,1:k-1) and L_(i,1:k-1) (i.e. swap row k with row i in L matrix, but only between columns 1 to k-1)
		L.ExchangeRows(k, i, 0, k-1);		

		//Interchange P_(k,:) and P (i,:) (i.e. swap row k with row i in P matrix, across all columns)
		P.ExchangeRows(k, i);
		
		//looping through the row indices of L and U
		for (int j = k+1; j < dim; ++j) {  
	
			int index_jk = GetIndex(j,k);	//to index jk entries  
			int index_kk = GetIndex(k,k);   //to index kk entries
			
			//to check if diagonal entries in U (aka pivots) are zero
			if (U.data[index_kk] == 0) {
				cerr << endl;
				cerr << "Error:  In Algorithm 2.2 - Gaussian Elimination with Partial Pivoting ..." << endl;
				cerr << "\tA zero pivot value has occurred." << endl;	
				cerr << "\tAlgorithm failed." << endl;
				cerr << endl;
				SquareMatrix error;   //use the default constructor to return the minimal version of SquareMatrix
				return error;
			}

			//if not then proceed with algorithm  (NB - no tmp data holders required here)
			else {
				//to alter L_jk value					
				L.data[index_jk] = U.data[index_jk]/U.data[index_kk];
			
				//to alter U values in row j from columns k to dim  (swapped i for a here to avoid ambiguity)
				for (int a = k; a < dim; ++a) {		
					int index_ja = GetIndex(j,a);	//to index j,k:m entries
					int index_ka = GetIndex(k,a);	//to index k,k:m entries			
					U.data[index_ja] -=  (L.data[index_jk]*U.data[index_ka]);
				}			
			}
		}
	}

	if (output == 'L') {
		return L;
	}
	if (output == 'U') {
		return U;
	}
	if (output == 'P') {
		return P;
	}
}

// Part 3. 

// 7. Forward Substitution - Algorithm 3.1

Matrix SquareMatrix::ForwardSub (const Matrix& b) const {
	
	//Check that invoking matrix is lower rectangular - error #1
	int dim = this->noOfCols;	
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			int index= GetIndex (i, j);	
			if (i<j) {
				double tolerance = 1e-14;
				if (abs(this->data[index]) > tolerance ) {		
					cerr << endl;
					cerr << "Error:  In Algorithm 3.1 - Forward Substitution..." << endl;
					cerr << "\tMatrix L is not lower rectangular." << endl;	
					cerr << endl;			
					Matrix error;   //use the default constructor to return the minimal version of Matrix
					return error;
				}
			}
		}
	}
	
	//Check that matrix b is of appropriate length - error #2
	int rows_b = b.GetNoOfRows();	
	if (this->noOfRows != rows_b) {
		cerr << endl;
		cerr << "Error:  In Algorithm 3.1 - Forward Substitution..." << endl;
		cerr << "\tNumber of rows in b does not match dimensions of L." << endl;
		cerr << endl;			
		Matrix error;   //use the default constructor to return the minimal version of Matrix
		return error;
	}

	//check that none of the diagonal entries in L are zero - error #3
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			int index= GetIndex (i, j);	
			if (i==j) {
				double tolerance = 1e-14;
				if (abs(this->data[index]) < tolerance ) {		
					cerr << endl;
					cerr << "Error:  In Algorithm 3.1 - Forward Substitution..." << endl;
					cerr << "\tMatrix L cannot have a zero valued diagonal entry." << endl;	
					cerr << endl;			
					Matrix error;   //use the default constructor to return the minimal version of Matrix
					return error;
				}
			}
		}
	}
	
	//once error checks are complete refer function to class Matrix for non-square matrix arithmetic and to access manipulate protected values in b
	return Matrix::ForwardSub(b);

}

// 8. Backward Substitution - Algorithm 3.2

Matrix SquareMatrix::BackwardSub (const Matrix& y) const {
	
	//Check that invoking matrix is upper rectangular - error #1
	int dim = this->noOfCols;	
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			int index= GetIndex (i, j);	
			if (i>j) {
				double tolerance = 1e-14;
				if (abs(this->data[index]) > tolerance ) {		
					cerr << endl;
					cerr << "Error:  In Algorithm 3.2 - Backward Substitution..." << endl;
					cerr << "\tMatrix U is not upper rectangular." << endl;	
					cerr << endl;			
					Matrix error;   //use the default constructor to return the minimal version of Matrix
					return error;
				}
			}
		}
	}
	
	//Check that matrix b is of appropriate length - error #2
	int rows_y = y.GetNoOfRows();	
	if (this->noOfRows != rows_y) {
		cerr << endl;
		cerr << "Error:  In Algorithm 3.2 - Backward Substitution..." << endl;
		cerr << "\tNumber of rows in y does not match dimensions of U." << endl;
		cerr << endl;			
		Matrix error;   //use the default constructor to return the minimal version of Matrix
		return error;
	}

	//check that none of the diagonal entries in U are zero - error #3
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			int index= GetIndex (i, j);	
			if (i==j) {
				double tolerance = 1e-14;
				if (abs(this->data[index]) < tolerance ) {		
					cerr << endl;
					cerr << "Error:  In Algorithm 3.2 - Backward Substitution..." << endl;
					cerr << "\tMatrix U cannot have a zero valued diagonal entry." << endl;	
					cerr << endl;			
					Matrix error;   //use the default constructor to return the minimal version of Matrix
					return error;
				}
			}
		}
	}	

	//once error checks are complete refer function to class Matrix for non-square matrix arithmetic and to access protected values in y
	return Matrix::BackwardSub(y);

}

// 9. Solving a system of linear equations in which the number of equations is equal to the number of unknowns, i.e. a uniquely determined solution exists. 

Matrix SquareMatrix::SolveUnique (Matrix& b) const {		
	
	//Check that matrix b is of appropriate length - error #1
	int rows_b = b.GetNoOfRows();	
	if (this->noOfRows != rows_b) {
		cerr << endl;
		cerr << "Error:  In function SolveUnique..." << endl;
		cerr << "\tNumber of rows in b does not match dimensions of A." << endl;
		cerr << endl;			
		Matrix error;   //use the default constructor to return the minimal version of Matrix
		return error;
	}

	//decompose A into L, U and P matrices
	SquareMatrix L = this->LUDecompositionTwo('L');
	//check that LU decomposition has occured succesfully - error #2
	if (L.noOfRows == 0) {
		//error message from LUDecompositionTwo will be displayed			
		Matrix error;  
		return error;
	}
	SquareMatrix U = this->LUDecompositionTwo('U');
	SquareMatrix P = this->LUDecompositionTwo('P');

	//pass these matrices to class Matrix for non-square matrix arithmetic and to access and manipulate protected values in b
	return Matrix::SolveUnique(L, U , P, b);
}
