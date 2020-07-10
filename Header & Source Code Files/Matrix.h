/*
 * Matrix.h
 *  
 *  Created on: 21.12.2017
 *      Author: Edward James - Medical Imaging CDT MRes, 2017-2018. 
 */

// Header file for class Matrix - provides support for arbitrarily sized m-by-n matrices

// Declares class data type called Matrix.

// Include header guards to avoid multiple attempts at defining class. 

#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>

using namespace std;

class Matrix {

// Part 1.

protected:
	int noOfRows; //stores the number of rows
	int noOfCols; //stores the number of columns 
	double* data; //stores the address to the 1-D array of the matrix entries arranged column-wise

	// getter	
	int GetIndex (const int& rowIdx, const int& columnIdx) const; // determines the position (index) along 'data' of a matrix entry, in the row and column specified by 'rowIdx' and 'columnIdx', respectively.

public:
	
// 1. ConstructorsTesting

	// constructors
	Matrix (); // default constructor -  must be defined as other constructors have been defined here
	Matrix (const int& noOfRows_rhs, const int& noOfCols_rhs); // standard custom constructor
	Matrix (const Matrix& input); // default copy constructor

	// named static constructors:
	static Matrix Zeros(const int& noOfRows_rhs, const int& noOfCols_rhs);
	static Matrix Ones (const int& noOfRows_rhs, const int& noOfCols_rhs);

	// default destructor 
	virtual ~Matrix ();

	// printing with ostream
	friend ostream& operator<< (ostream& out, const Matrix& rhs);

// 2. AssignmentTesting
	
	//default assignment operator
	Matrix& operator= (const Matrix& rhs);

// 3. ToeplitzTesting & ToeplitzTestingHelper
	
	static Matrix Print (const double* const data_rhs, const int& noOfRows_rhs, const int& noOfCols_rhs);
	static Matrix Toeplitz (const double* const column, const int& noOfRows_rhs, const double* const row, const int& noOfCols_rhs); //creates a new Toeplitz matrix

// 4. TransposeTesting

	static Matrix Transpose (const Matrix& input);  //static function - creates a new transpose matrix
	void Transpose ();  //non-static function - converts matrix into tranpose of itself

// 5. MultiplicationTesting

	friend Matrix operator* (const Matrix& lhs, const Matrix& rhs); //non-member multiplication function - creates new matrix which is result of lhs * rhs
	Matrix& operator*= (const Matrix& rhs); // member self-modifying multiplication function - converts input (lhs) matrix into result of lhs * rhs

// 6. RowColumnExchangeTesting
	
	Matrix& ExchangeRows (const int& row1, const int& row2); //swaps row1 and row2, for all columns in matrix
	Matrix& ExchangeRows (const int& row1, const int& row2, const int& col1, const int& col2); //swaps row1 and row2, but only between col1 and col2
	Matrix& ExchangeColumns (const int& col1, const int& col2); //swaps col1 and col2, for all rows in matrix
	Matrix& ExchangeColumns (const int& col1, const int& col2, const int& row1, const int& row2); //swaps col1 and col2, but only between row1 and row2	


// 7. OtherTesting

	void Zeros (); // (sets every entry of matrix to zero)
	void Ones (); // (sets every entry of matrix to one)
	int GetNoOfRows () const; // (find out number of rows in matrix)
	int GetNoOfColumns () const; // (find out number of columns in matrix)
	double GetEntry (const int& rowIdx, const int& columnIdx) const; // (find out value of particular entry (i,j) in Matrix, starting indexing at (0,0) in upper left corner of matrix)

// Part 2.

	static Matrix Test (const double* data_rhs, const int& noOfRows_rhs, const int& noOfCols_rhs); //for constructing specified test matrices for testing in Parts 2 & 3	
	static Matrix PrintMATLAB (const Matrix& rhs); //Specifically to print matrices in default MATLAB style to fixed 4.d.p precision for Parts 2 & 3
	
// Part 3. 

// 8. Forward Substitution - Algorithm 3.1

	virtual Matrix ForwardSub (const Matrix& b) const; 
		//Invoking matrix L is a lower triangular m-by-m square matrix (from SquareMatrix class function). 
		//Input b is a m-by-1 column vector. 
		//Returns a m-by-1 column vector y, such that Ly=b. 

// 9. Backward Substitution - Algorithm 3.2

	virtual Matrix BackwardSub (const Matrix& y) const; 
		//Invoking matrix U is an upper triangular m-by-m square matrix (from SquareMatrix class function). 
		//Input y is a m-by-1 column vector. 
		//Returns a m-by-1 column vector x, such that Ux=y. 

// 10. Solving a system of linear equations in which the number of equations is equal to the number of unknowns, i.e. a uniquely determined solution exists. 

	Matrix SolveUnique (const Matrix& L, const Matrix& U, const Matrix& P, Matrix& b) const;
		//Function is invoked by a m-by-m square matrix A (from SquareMatrix class function), which represents the matrix of coefficients
		//Input matrices L, U and P are produced in class SquareMatrix from LU decomposition of A using Algorithm 2.2. All 3 are m-by-m square matrices. 
		//Input b is a m-by-1 column vector.
		//Returns a m-by-1 column vector x, such that Ax=b.

// 11. Solving an overdetermined system of linear equations in which the number of equations is greater than the number of unknowns, i.e. a least squares solution is produced. 

	Matrix LeastSquaresHelper (const Matrix& b, const string& output) const;
		//Function is invoked by a non-square m-by-n matrix A, which represents the matrix of coefficients 
		//Input b is a m-by-1 column vector.
		//Returns n-by-n square matrix A'*A if output = 'A'*A', or n-by-1 column vector A'*b if output = 'A'*b'. 

};

#endif /* MATRIX_H_ */
