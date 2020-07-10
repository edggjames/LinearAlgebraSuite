/*
 * SquareMatrix.h
 *  
 *  Created on: 21.12.2017
 *      Author: Edward James - Medical Imaging CDT MRes, 2017-2018. 
 */

// Header file for class SquareMatrix - extends class Matrix to provide specialised support for square m-by-m matrices.

// Declares class data type called SquareMatrix.

// Include header guards to avoid multiple attempts at defining class.  

#ifndef SQUAREMATRIX_H_
#define SQUAREMATRIX_H_

#include "Matrix.h"

class SquareMatrix : public Matrix { // class SquareMatrix inherits publicly from class Matrix

private:

	void checkSquare () const;  //to check if user is attempting to make square copy of unsquare matrix

public:

// Part 2. 

// 1. ConstructorsTesting

	//constructors
	SquareMatrix (); //default constructor  
	SquareMatrix (const int& dim); // standard custom constructor
	SquareMatrix (const SquareMatrix& input); // default copy constructor #1  //returns SquareMatrix copy of SquareMatrix
	SquareMatrix (const Matrix& input); // default copy constructor #2        //returns SquareMatrix copy of Matrix

	// named static constructors:
	static SquareMatrix Zeros (const int& dim); 
	static SquareMatrix Ones (const int& dim); 
	static SquareMatrix Eye (const int& dim);
	static SquareMatrix Test (const double* data_rhs, const int& dim); //for constructing specified test matrices for testing

	// default destructor
	~SquareMatrix ();

// 2. ToeplitzTesting & ToeplitzTestingHelper

	static SquareMatrix Toeplitz (const double* const row, const int& dim); //creates a new Toeplitz matrix when only first row is specified
	static SquareMatrix Toeplitz (const double* const column, const double* const row, const int& dim); //creates a new Toeplitz matrix when both first column and first row are specified

// 3. TransposeTesting

	//no new function declarations required.

// 4. TriangularExtractionTesting
	
	//extract upper triangular part of square matrix
	SquareMatrix TriUpper () const;  // non-static member function

	//extract lower triangular part of square matrix
	SquareMatrix TriLower () const;  // non-static member function

// 5. LUDecompositionTesting - Algorithm 2.1 - Gaussian Elimination without Pivoting

	SquareMatrix LUDecompositionOne (const char& output) const;
		//Takes a non-singular square matrix A (i.e. with non-zero determinant) as input (i.e. invoking square matrix).
		//Returns lower triangular matrix L if output = 'L', or upper triangular matrix U if output = 'U', such that A=LU.	

// 6. LUDecompositionTesting - Algorithm 2.2 - Gaussian Elimination with Partial Pivoting

	SquareMatrix LUDecompositionTwo (const char& output) const;
		//Takes a non-singular square matrix A (i.e. with non-zero determinant) as input (i.e. invoking square matrix).
		//Returns lower triangular matrix L if output = 'L', or upper triangular matrix U if output = 'U', or permutation matrix P if output = 'P', such that PA=LU.

// Part 3. 

// 7. Forward Substitution - Algorithm 3.1

	Matrix ForwardSub (const Matrix& b) const;
		//Invoking matrix L is a lower triangular m-by-m square matrix. 
		//Input b is a m-by-1 column vector. 
		//Returns a m-by-1 column vector y, such that Ly=b (via a call to the Matrix class function). 

// 8. Backward Substitution - Algorithm 3.2

	Matrix BackwardSub (const Matrix& y) const; 
		//Invoking matrix U is an upper triangular m-by-m square matrix. 
		//Input y is a m-by-1 column vector. 
		//Returns a m-by-1 column vector x, such that Ux=y (via a call to the Matrix class function).

// 9. Solving a system of linear equations in which the number of equations is equal to the number of unknowns, i.e. a uniquely determined solution exists. 

	Matrix SolveUnique (Matrix& b) const;
		//Function is invoked by a m-by-m square matrix A, which represents the matrix of coefficients
		//Input b is a m-by-1 column vector.
		//Returns a m-by-1 column vector x, such that Ax=b (via a call to the Matrix class function).

};

#endif /* SQUAREMATRIX_H_ */
