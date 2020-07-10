	/*
 * SquareMatrixTest.cpp
 *
 *  Created on: 21.12.2017
 *      Adapted from code provided by: Gary Hui Zhang (gary.zhang@ucl.ac.uk)
 */

// A program to demonstrate and test the Matrix class and the SquareMatrix class

#include "Matrix.h"
#include "SquareMatrix.h"
#include <iostream>
#include <cstdlib>  //need this for 'system' functions


using namespace std;

void ConstructorsTesting () {
	cout << "Testing the SquareMatrix constructors:" << endl;
	cout << endl;

	cout << "Case 1: Creating a square matrix of zeros with the standard constructor:" << endl;
	{
		SquareMatrix matrix(4);  //4 is dim in example.m
		cout << matrix << endl;		
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: Creating a square matrix of zeros with the static Zeros function: " << endl;
	{
		SquareMatrix matrix = SquareMatrix::Zeros(4);
		cout << matrix << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}
	
	cout << "Case 3: Creating a square matrix of ones with the static Ones function: " << endl;
	{
		SquareMatrix matrix1 = SquareMatrix::Ones(4);
		cout << matrix1 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 4: Creating an identity matrix with the static Eye function:" << endl;
	{
		SquareMatrix matrix2 = SquareMatrix::Eye(4);
		cout << matrix2 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 5: Copying an identity matrix with a copy constructor:" << endl;
	{
		SquareMatrix matrix2 = SquareMatrix::Eye(4);
		cout << endl;
		cout << "The input matrix = " << endl;
		cout << matrix2 << endl;
		SquareMatrix copy_matrix2 (matrix2); 
		cout << "The copy = " << endl;
		cout << copy_matrix2 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 6: Constructing a user defined matrix with the test constructor:" << endl;
	{
		double vector[16] = {2, 4, 8, 6, 1, 3, 7, 7, 1, 3, 9, 9, 0, 1, 5, 8};  
		SquareMatrix matrix = SquareMatrix::Test(vector, 4);
		cout << matrix << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}
	
	return;
}

//'first row only' overloaded variety
void ToeplitzTestingHelper (const double* const row, const int& dim, const double* const expected) {
	cout << "The 1st row = " << endl;
	Matrix::Print(row, 1, dim);
	cout << endl;
	cout << "The matrix created by the toeplitz function in MATLAB = " << endl;
	Matrix::Print(expected, dim, dim);
	cout << endl;
	SquareMatrix matrix3 = SquareMatrix::Toeplitz(row, dim);
	cout << "The matrix created by SquareMatrix::Toeplitz = " << endl;
	cout << matrix3 << endl;

	return;
}

//'first column and first row' overloaded variety
void ToeplitzTestingHelper (const double* const column, const double* const row, const int& dim, const double* const expected) {
	cout << "The 1st column = " << endl;
	Matrix::Print(column, dim, 1);
	cout << endl;	
	cout << "The 1st row = " << endl;
	Matrix::Print(row, 1, dim);
	cout << endl;
	cout << "The matrix created by the toeplitz function in MATLAB = " << endl;
	Matrix::Print(expected, dim, dim);
	cout << endl;
	SquareMatrix matrix4 = SquareMatrix::Toeplitz(column, row, dim);
	cout << "The matrix created by SquareMatrix::Toeplitz = " << endl;
	cout << matrix4 << endl;

	return;
}

void ToeplitzTesting () {
	cout << "Testing the static functions SquareMatrix::Toeplitz:" << endl;
	cout << endl;

	cout << "Case 1: When only the first row is specified:" << endl;
	cout << endl;

	{
		double row[4] = {4, 3, 2, 1};  //aka vec1
		// matrix should be stored, in 1-D, column-wise
		//  4     3     2     1
		//  3     4     3     2
		//  2	  3	4     3
		//  1 	  2     3     4
		double expected[16] = {4, 3, 2, 1, 3, 4, 3, 2, 2, 3, 4, 3, 1, 2, 3, 4};
		ToeplitzTestingHelper(row, 4, expected);
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: When both the first column and the first row are specified, and the" << endl;
	cout << "        first elements are equal:" << endl;
	cout << endl;

	{
		double column[4] = {2, 1, 0, -1};
		double row[4] = {2, 3, 4, 5};
		// matrix should be stored, in 1-D, column-wise
		//  2     3     4     5
		//  1     2     3     4
		//  0     1     2     3
		// -1     0     1     2
		double expected[16] = {2, 1, 0, -1, 3, 2, 1, 0, 4, 3, 2, 1, 5, 4, 3, 2};
		ToeplitzTestingHelper(column, row, 4, expected);
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: When both the first column and the first row are specified, and the" << endl;
	cout << "        first elements of each are not equal:" << endl;
	cout << endl;

	{
		double column[4] = {1, 2, 3, 4}; //aka vec2		
		double row[4] = {4, 3, 2, 1};	 //aka vec1	
		// matrix should be stored, in 1-D, column-wise
		//  1     3     2     1
		//  2     1     3     2
		//  3     2     1     3	
		//  4	  3 	2     1
		double expected[16] = {1, 2, 3, 4, 3, 1, 2, 3, 2, 3, 1, 2, 1, 2, 3, 1};
		ToeplitzTestingHelper(column, row, 4, expected);
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	return;
}

void TransposeTesting () {
	cout << "Testing the class Matrix Transpose functions for SquareMatrix objects:" << endl;
	cout << endl;	

	cout << "Case 1: the non-static self-modifying Transpose function:" << endl;
	cout << endl;	

	{
		// the same matrix as in ToeplitzTesting case 3
		double column[4] = {1, 2, 3, 4};  //aka vec2		
		double row[4] = {4, 3, 2, 1};     //aka vec1
		SquareMatrix matrix4 = SquareMatrix::Toeplitz(column, row, 4);
		cout << "The original Matrix = " << endl;
		cout << matrix4 << endl;
		matrix4.Transpose();
		SquareMatrix matrix5 = matrix4;
		cout << "The transposed version = " << endl;
		cout << matrix5 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: the static Transpose function:" << endl;
	cout << endl;

	{
		// the same matrix as in ToeplitzTesting case 3
		double column[4] = {1, 2, 3, 4}; //aka vec2		
		double row[4] = {4, 3, 2, 1};	 //aka vec1
		SquareMatrix matrix4 = SquareMatrix::Toeplitz(column, row, 4);
		cout << "The original Matrix = " << endl;
		cout << matrix4 << endl;
		SquareMatrix matrix5 = Matrix::Transpose(matrix4);
		cout << "The transposed version = " << endl;		
		cout << matrix5 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}
	
	return;	

}

void TriangularExtractionTesting () { 
	cout << "Testing the triangular extraction functions of class SquareMatrix:" << endl;
	cout << endl;

	cout << "Case 1: the upper triangular extraction non-static function:" << endl;

	{
		// the same matrix as in ToeplitzTesting case 3
		double column[4] = {1, 2, 3, 4}; //aka vec2		
		double row[4] = {4, 3, 2, 1};	 //aka vec1
		cout << endl;		
		SquareMatrix matrix4 = SquareMatrix::Toeplitz(column, row, 4);
		cout << "The original square matrix = " << endl;
		cout << matrix4 << endl;
		SquareMatrix matrix6 = matrix4.TriUpper();
		cout << "The square matrix after upper triangular extraction = " << endl;
		cout << matrix6 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}
	
	cout << "Case 2: the lower triangular extraction non-static function:" << endl;

	{
		// the same matrix as in ToeplitzTesting case 3
		double column[4] = {1, 2, 3, 4}; //aka vec2		
		double row[4] = {4, 3, 2, 1};	 //aka vec1
		cout << endl;		
		SquareMatrix matrix4 = SquareMatrix::Toeplitz(column, row, 4);
		cout << "The original square matrix = " << endl;
		cout << matrix4 << endl;
		SquareMatrix matrix7 = matrix4.TriLower();
		cout << "The square matrix after lower triangular extraction = " << endl;
		cout << matrix7 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	return;

}

void LUDecompositionTestingOne() {
	cout << "Testing the LU decomposition of a square matrix using Algorithm 2.1, Gaussian" << endl;
	cout << "elimination without pivoting:" << endl;
	cout << endl;

	cout << "Case 1: Matrix A is non-singular and all entries on main diagonal are non-zero:" << endl;
	cout << endl;

	{
		// the same matrix as in ToeplitzTesting case 1
		double row[4] = {4, 3, 2, 1};  //aka vec1
		SquareMatrix matrix3 = SquareMatrix::Toeplitz(row, 4);
		cout << "The input Matrix (A) = " << endl;
		cout << matrix3 << endl;
		SquareMatrix lmatrix1 = matrix3.LUDecompositionOne('L');
		cout << "The L Matrix = " << endl;
		Matrix::PrintMATLAB(lmatrix1);   
 		SquareMatrix umatrix1 = matrix3.LUDecompositionOne('U');
		cout << "The U Matrix = " << endl;
		Matrix::PrintMATLAB(umatrix1);   
		cout << "Checking that (L*U = A):" << endl;
		cout << lmatrix1 * umatrix1 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: Matrix A is non-singular and 2nd, 3rd and 4th entries on main diagonal" << endl;
	cout << "        are all zero:" << endl;
	cout << endl;

	{
		double data[16] = {6, 3, 6, 7, 9, 0, 3, 8, 7, 1, 0, 9, 4, 2, 1, 0};  
		SquareMatrix matrix = SquareMatrix::Test(data, 4);
		cout << "The input Matrix (A) = " << endl;
		cout << matrix << endl;
		SquareMatrix lmatrix = matrix.LUDecompositionOne('L');
		cout << "The L Matrix = " << endl;
		Matrix::PrintMATLAB(lmatrix);   
 		SquareMatrix umatrix = matrix.LUDecompositionOne('U');
		cout << "The U Matrix = " << endl;
		Matrix::PrintMATLAB(umatrix);   
		cout << "Checking that (L*U = A):" << endl;
		cout << lmatrix * umatrix << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: Matrix A is non-singular and only 1st entry on main diagonal is zero:" << endl;
	cout << endl;

	{
		double data[16] = {0, 3, 6, 7, 9, 2, 3, 8, 7, 1, 2, 9, 4, 2, 1, 5};  
		SquareMatrix matrix = SquareMatrix::Test(data, 4);
		cout << "The input Matrix (A) = " << endl;
		cout << matrix << endl;
		cout << "Attempting to decompose into L and U matrices ... " << endl;		
		SquareMatrix lmatrix1 = matrix.LUDecompositionOne('L');
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	return;	

}

void LUDecompositionTestingTwo() {
	cout << "Testing the LU decomposition of a square matrix using Algorithm 2.2, Gaussian" << endl;
	cout << "elimination with partial pivoting:" << endl;
	cout << endl;

	cout << "Case 1: Matrix A is non-singular and magnitude of each entry on main diagonal is" << endl;
	cout << "        maximal within its column. Partial pivoting is not required: " << endl;
	cout << endl;

	{
		// the same matrix as in ToeplitzTesting case 1
		double row[4] = {4, 3, 2, 1};  //aka vec1
		SquareMatrix matrix3 = SquareMatrix::Toeplitz(row, 4);
		cout << "The input Matrix (A) = " << endl;
		cout << matrix3 << endl;
		SquareMatrix lmatrix1 = matrix3.LUDecompositionTwo('L');
		cout << "The L Matrix = " << endl;
		Matrix::PrintMATLAB(lmatrix1);   
 		SquareMatrix umatrix1 = matrix3.LUDecompositionTwo('U');
		cout << "The U Matrix = " << endl;
		Matrix::PrintMATLAB(umatrix1);   
		SquareMatrix pmatrix1 = matrix3.LUDecompositionTwo('P');
		cout << "The P Matrix = " << endl;
		cout << pmatrix1 << endl;   
		cout << "Checking that (L*U = P*A): " << endl;
		cout << endl;
		cout << "L*U = " << endl;
		cout << lmatrix1 * umatrix1 << endl;
		cout << "P*A = " << endl;		
		cout << pmatrix1 * matrix3 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}
	
	cout << "Case 2: Matrix A is non-singular and magnitude of each entry on main diagonal is" << endl;
	cout << "        not maximal within its column. Partial pivoting is required: " << endl;
	cout << endl;

	{
		// the same matrix as in ToeplitzTesting case 3
		double column[4] = {1, 2, 3, 4}; //aka vec2		
		double row[4] = {4, 3, 2, 1};	 //aka vec1		
		SquareMatrix matrix4 = SquareMatrix::Toeplitz(column, row, 4);
		cout << "The input Matrix (A) = " << endl;
		cout << matrix4 << endl;
		SquareMatrix lmatrix2 = matrix4.LUDecompositionTwo('L');
		cout << "The L Matrix = " << endl;
		Matrix::PrintMATLAB(lmatrix2);   
 		SquareMatrix umatrix2 = matrix4.LUDecompositionTwo('U');
		cout << "The U Matrix = " << endl;
		Matrix::PrintMATLAB(umatrix2);   
		SquareMatrix pmatrix2 = matrix4.LUDecompositionTwo('P');
		cout << "The P Matrix = " << endl;
		cout << pmatrix2 << endl;   
		cout << "Checking that (L*U = P*A): " << endl;
		cout << endl;
		cout << "L*U = " << endl;
		cout << lmatrix2 * umatrix2 << endl;
		cout << "P*A = " << endl;		
		cout << pmatrix2 * matrix4 << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: Matrix A is non-singular and all diagonal entries are zero:" << endl;
	cout << endl;

	{
		double data[16] = {0, 3, 6, 7, 9, 0, 3, 8, 7, 1, 0, 9, 4, 2, 1, 0};  
		SquareMatrix matrix = SquareMatrix::Test(data, 4);;
		cout << "The input Matrix (A) = " << endl;
		cout << matrix << endl;
		SquareMatrix lmatrix = matrix.LUDecompositionTwo('L');
		cout << "The L Matrix = " << endl;
		Matrix::PrintMATLAB(lmatrix);   
 		SquareMatrix umatrix = matrix.LUDecompositionTwo('U');
		cout << "The U Matrix = " << endl;
		Matrix::PrintMATLAB(umatrix);   
		SquareMatrix pmatrix = matrix.LUDecompositionTwo('P');
		cout << "The P Matrix = " << endl;
		cout << pmatrix << endl;   
		cout << "Checking that (L*U = P*A): " << endl;
		cout << endl;
		cout << "L*U = " << endl;
		cout << lmatrix * umatrix << endl;
		cout << "P*A = " << endl;		
		cout << pmatrix * matrix << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 4: Matrix A is non-singular, all diagonal entries are zero, and first three" << endl;
	cout << "        entries of first column are zero:" << endl;
	cout << endl;

	{
		double data[16] = {0, 0, 0, 7, 9, 0, 3, 8, 7, 1, 0, 9, 4, 2, 1, 0};  
		SquareMatrix matrix = SquareMatrix::Test(data, 4);;
		cout << "The input Matrix (A) = " << endl;
		cout << matrix << endl;
		SquareMatrix lmatrix = matrix.LUDecompositionTwo('L');
		cout << "The L Matrix = " << endl;
		Matrix::PrintMATLAB(lmatrix);   
 		SquareMatrix umatrix = matrix.LUDecompositionTwo('U');
		cout << "The U Matrix = " << endl;
		Matrix::PrintMATLAB(umatrix);   
		SquareMatrix pmatrix = matrix.LUDecompositionTwo('P');
		cout << "The P Matrix = " << endl;
		cout << pmatrix << endl;   
		cout << "Checking that (L*U = P*A): " << endl;
		cout << endl;
		cout << "L*U = " << endl;
		cout << lmatrix * umatrix << endl;
		cout << "P*A = " << endl;		
		cout << pmatrix * matrix << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 5: Matrix A is singular, all diagonal entries are zero, and all entries" << endl;
	cout << "        of first column are zero:" << endl;
	cout << endl;

	{
		double data[16] = {0, 0, 0, 0, 9, 0, 3, 8, 7, 1, 0, 9, 4, 2, 1, 0};  
		SquareMatrix matrix = SquareMatrix::Test(data, 4);
		cout << "The input Matrix (A) = " << endl;
		cout << matrix << endl;
		cout << "Attempting to decompose into L, U and P matrices ... " << endl;	
 		SquareMatrix umatrix = matrix.LUDecompositionTwo('L');
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	return;
}

void ForwardSubstitution () { 

	cout << "Testing Algorithm 3.1 - Forward Substitution. To solve for y in L*y = b:" << endl;
	cout << endl;

	cout << "Case 1: Successful Forward Substitution: " << endl;
	cout << endl;

	{
		double matrix[9] = {1, 3, 4, 0, -1, 1, 0, 0, -3};  
		SquareMatrix L = SquareMatrix::Test(matrix, 3);
		cout << "The lower triangular matrix, L = " << endl;
		cout << L << endl;
		double vector[3] = {16, 43, 57};  
		Matrix b = Matrix::Test(vector, 3, 1);
		cout << "The column vector, b  = " << endl;
		cout << b << endl;
		cout << "The result of the forward substitution, the column vector y = " << endl;		
		Matrix y = L.ForwardSub(b);
		cout << y << endl;
		cout << "Checking that (L*y = b): " << endl;
		cout << L * y << endl;   
 		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: Testing error message #1:" << endl;
	cout << endl;

	{
		double matrix[9] = {1, 3, 4, 0, -1, 1, 0, 1, -3};  
		SquareMatrix L = SquareMatrix::Test(matrix, 3);
		cout << "The input matrix, 'L' = " << endl;
		cout << L << endl;
		double vector[3] = {16, 43, 57};  
		Matrix b = Matrix::Test(vector, 3, 1);
		cout << "The column vector, b  = " << endl;
		cout << b << endl;
		cout << "Attempting forward substitution ... " << endl;
		Matrix y = L.ForwardSub(b); 
 		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: Testing error message #2:" << endl;
	cout << endl;

	{
		double matrix[9] = {1, 3, 4, 0, -1, 1, 0, 0, -3};  
		SquareMatrix L = SquareMatrix::Test(matrix, 3);
		cout << "The lower triangular matrix, L = " << endl;
		cout << L << endl;
		double vector[4] = {16, 43, 57, 29};  
		Matrix b = Matrix::Test(vector, 4, 1);
		cout << "The column vector, b  = " << endl;
		cout << b << endl;
		cout << "Attempting forward substitution ... " << endl;
		Matrix y = L.ForwardSub(b); 
 		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 4: Testing error message #3:" << endl;
	cout << endl;

	{
		double matrix[9] = {1, 3, 4, 0, -1, 1, 0, 0, 0};  
		SquareMatrix L = SquareMatrix::Test(matrix, 3);
		cout << "The lower triangular matrix, L = " << endl;
		cout << L << endl;
		double vector[3] = {16, 43, 57};  
		Matrix b = Matrix::Test(vector, 3, 1);
		cout << "The column vector, b  = " << endl;
		cout << b << endl;
		cout << "Attempting forward substitution ... " << endl;
		Matrix y = L.ForwardSub(b);   
 		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}
	
	return;	

}

void BackwardSubstitution () { 

	cout << "Testing Algorithm 3.2 - Backward Substitution. To solve for x in U*x = y:" << endl;
	cout << endl;

	cout << "Case 1: Successful Backward Substitution: " << endl;
	cout << endl;

	{
		double matrix[9] = {1, 0, 0, 2, 1, 0, 3, 1, 1};  
		SquareMatrix U = SquareMatrix::Test(matrix, 3);
		cout << "The upper triangular matrix, U = " << endl;
		cout << U << endl;
		double vector[3] = {16, 5, 4};  
		Matrix y = Matrix::Test(vector, 3, 1);
		cout << "The column vector, y  = " << endl;
		cout << y << endl;
		cout << "The result of the backward substitution, the column vector x = " << endl;
		Matrix x = U.BackwardSub(y);
		cout << x << endl;
		cout << "Checking that (U*x = y): " << endl;
		cout << U * x << endl;   
 		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: Testing error message #1:" << endl;
	cout << endl;

	{
		double matrix[9] = {1, 0, 1, 2, 1, 0, 3, 1, 1}; 
		SquareMatrix U = SquareMatrix::Test(matrix, 3);
		cout << "The input matrix, 'U' = " << endl;
		cout << U << endl;
		double vector[3] = {16, 43, 57};  
		Matrix y = Matrix::Test(vector, 3, 1);
		cout << "The column vector, y  = " << endl;
		cout << y << endl;
		cout << "Attempting backward substitution ... " << endl;
		Matrix x = U.BackwardSub(y); 
 		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: Testing error message #2:" << endl;
	cout << endl;

	{
		double matrix[9] = {1, 0, 0, 2, 1, 0, 3, 1, 1};  
		SquareMatrix U = SquareMatrix::Test(matrix, 3);
		cout << "The upper triangular matrix, U = " << endl;
		cout << U << endl;
		double vector[4] = {16, 43, 57, 29};  
		Matrix y = Matrix::Test(vector, 4, 1);
		cout << "The column vector, y  = " << endl;
		cout << y << endl;
		cout << "Attempting backward substitution ... " << endl;
		Matrix x = U.BackwardSub(y); 
 		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 4: Testing error message #3:" << endl;
	cout << endl;

	{
		double matrix[9] = {0, 0, 0, 2, 1, 0, 3, 1, 1};  
		SquareMatrix U = SquareMatrix::Test(matrix, 3);
		cout << "The upper triangular matrix, U = " << endl;
		cout << U << endl;
		double vector[3] = {16, 43, 57};  
		Matrix y = Matrix::Test(vector, 3, 1);
		cout << "The column vector, y  = " << endl;
		cout << y << endl;
		cout << "Attempting backward substitution ... " << endl;
		Matrix x = U.BackwardSub(y); 
 		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	return;	

}

void UniqueSolution() {

	cout << "Solving a system of linear equations where a unique solution exists:" << endl;
	cout << endl;

	cout << "Case 1: Successful solution for a 4 x 4 matrix:" << endl;
	cout << endl;

	{
		double matrix[16] = {2, 4, 8, 6, 1, 3, 7, 7, 1, 3, 9, 9, 0, 1, 5, 8};  
		SquareMatrix A = SquareMatrix::Test(matrix, 4);
		cout << "The matrix of coefficients, A = " << endl;
		cout << A << endl;
		double vector[4] = {1, 0, 0, 0};  
		Matrix b = Matrix::Test(vector, 4, 1);
		cout << "The column vector, b  = " << endl;
		cout << b << endl;
		Matrix x = A.SolveUnique(b);	
		cout << "The result of solving for x in (A*x = b), the column vector x  = " << endl;				
		Matrix::PrintMATLAB(x); 
		cout << "Checking that (A*x = b): " << endl;
		cout << A * x << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}
	
	cout << "Case 2: Successful solution for a 5 x 5 matrix:" << endl;
	cout << endl;

	{
		double matrix[25] = {2, 4, 8, 6, 1, 3, 7, 7, 1, 3, 9, 9, 0, 1, 5, 8, 4, 6, 8, 2, 0, 9, 7, 3, 4};  
		SquareMatrix A = SquareMatrix::Test(matrix, 5);
		cout << "The matrix of coefficients, A = " << endl;
		cout << A << endl;
		double vector[5] = {6, 5, 3, 1, 4};
		Matrix b = Matrix::Test(vector, 5, 1);
		cout << "The column vector, b  = " << endl;
		cout << b << endl;
		Matrix x = A.SolveUnique(b);	
		cout << "The result of solving for x in (A*x = b), the column vector x  = " << endl;				
		Matrix::PrintMATLAB(x); 
		cout << "Checking that (A*x = b): " << endl;
		cout << A * x << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: Error message testing:" << endl;
	cout << endl;

	{
		double matrix[16] = {2, 4, 8, 6, 1, 3, 7, 7, 1, 3, 9, 9, 0, 1, 5, 8};  
		SquareMatrix A = SquareMatrix::Test(matrix, 4);
		cout << "The matrix of coefficients, A = " << endl;
		cout << A << endl;
		double vector[5] = {1, 0, 0, 0, 2};  
		Matrix b = Matrix::Test(vector, 5, 1);
		cout << "The column vector, b  = " << endl;
		cout << b << endl;
		cout << "Attempting to solve for x in (A*x = b) ... " << endl;
		Matrix x = A.SolveUnique(b);		
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 4: Error message testing:" << endl;
	cout << endl;

	{
		double data[16] = {0, 0, 0, 0, 9, 0, 3, 8, 7, 1, 0, 9, 4, 2, 1, 0};  
		SquareMatrix A = SquareMatrix::Test(data, 4);
		cout << "When the matrix of coefficients is singular. A = " << endl;
		cout << A << endl;
		double vector[4] = {1, 0, 0, 0};  
		Matrix b = Matrix::Test(vector, 4, 1);
		cout << "The column vector, b  = " << endl;
		cout << b << endl;
		cout << "Attempting to solve for x in (A*x = b) ... " << endl;
		Matrix x = A.SolveUnique(b);		
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	return;
}

void LeastSquareSolution() {

	cout << "Solving an overdetermined system of linear equations for a least squares" << endl;
	cout << "solution:" << endl;
	cout << endl;

	cout << "Case 1: Successful solution for a 4 x 2 matrix:" << endl;
	cout << endl;

	{
		double matrix[8] = {4, 3, 2, 1, 3, 4, 3, 2};  
		Matrix A = Matrix::Test(matrix, 4, 2);
		cout << "The m-by-n non-square matrix of coefficients, A = " << endl;
		cout << A << endl;
		double vector1[4] = {1, 2, 3, 4};  
		Matrix b = Matrix::Test(vector1, 4, 1);
		cout << "The m-by-1 column vector, b  = " << endl;
		cout << b << endl;
		SquareMatrix ATA = A.LeastSquaresHelper(b,"A'*A");
		cout << "Converting A into A'*A, a n-by-n matrix = " << endl;		
		Matrix ATb = A.LeastSquaresHelper(b, "A'*b");		
		cout << ATA << endl;
		cout << "Converting b into A'*b, a n-by-1 column vector = " << endl;
		cout << ATb << endl;
		Matrix x = ATA.SolveUnique(ATb);	
		cout << "The result of solving for x in (A'A*x = A'*b), the n-by-1 column vector x  = " << endl;		
		Matrix::PrintMATLAB(x); 
		cout << "Checking against the solution given in MATLAB for x = A\\b : " << endl;
		double vector2[2] = {-1.1724, 1.7241};
		Matrix expected = Matrix::Test(vector2, 2, 1);	
		Matrix::PrintMATLAB(expected);	
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: Successful solution for a 5 x 4 matrix:" << endl;
	cout << endl;

	{
		double matrix[20] = {4, 3, 2, 1, 3, 4, 3, 2, 2, 3, 4, 3, 1, 2, 3, 4, 4, 3, 2, 1};  
		Matrix A = Matrix::Test(matrix, 5, 4);
		cout << "The m-by-n non-square matrix of coefficients, A = " << endl;
		cout << A << endl;
		double vector1[5] = {1, 2, 3, 4, 5};  
		Matrix b = Matrix::Test(vector1, 5, 1);
		cout << "The m-by-1 column vector, b  = " << endl;
		cout << b << endl;
		SquareMatrix ATA = A.LeastSquaresHelper(b,"A'*A");
		cout << "Converting A into A'*A, a n-by-n matrix = " << endl;		
		Matrix ATb = A.LeastSquaresHelper(b, "A'*b");		
		cout << ATA << endl;
		cout << "Converting b into A'*b, a n-by-1 column vector = " << endl;
		cout << ATb << endl;
		Matrix x = ATA.SolveUnique(ATb);	
		cout << "The result of solving for x in (A'A*x = A'*b), the n-by-1 column vector x  = " << endl;		
		Matrix::PrintMATLAB(x); 
		cout << "Checking against the solution given in MATLAB for x = A\\b : " << endl;
		double vector2[4] = {-2.7453, 7.5652, -2.9689, -1.2236};
		Matrix expected = Matrix::Test(vector2, 4, 1);	
		Matrix::PrintMATLAB(expected);	
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: Testing warning message:" << endl;
	cout << endl;

	{

		double column[4] = {1, 2, 3, 4};		
		SquareMatrix A = SquareMatrix::Toeplitz(column,4);
		cout << "When the matrix of coefficients is square. A = " << endl;
		cout << A << endl;
		double vector1[4] = {1, 2, 3, 4};  
		Matrix b = Matrix::Test(vector1, 4, 1);
		cout << "The m-by-1 column vector, b  = " << endl;
		cout << b << endl;
		SquareMatrix ATA = A.LeastSquaresHelper(b,"A'*A");
		cout << "Converting A into A'*A, a n-by-n matrix = " << endl;		
		Matrix ATb = A.LeastSquaresHelper(b, "A'*b");		
		cout << ATA << endl;
		cout << "Converting b into A'*b, a n-by-1 column vector = " << endl;
		cout << ATb << endl;
		Matrix x = ATA.SolveUnique(ATb);	
		cout << "The result of solving for x in (A'A*x = A'*b), the n-by-1 column vector x  = " << endl;		
		Matrix::PrintMATLAB(x); 
		cout << "Checking that (A*x = b): " << endl;
		cout << A * x << endl;
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 4: Testing error message #1:" << endl;
	cout << endl;

	{
		double matrix[8] = {1, 1, 1, 1, 1, 2, 3, 4};  
		Matrix A = Matrix::Test(matrix, 2, 4);
		cout << "The m-by-n non-square matrix of coefficients, A = " << endl;
		cout << A << endl;
		double vector1[2] = {6, 5};  
		Matrix b = Matrix::Test(vector1, 2, 1);
		cout << "The m-by-1 column vector, b  = " << endl;
		cout << b << endl;
		cout << "Attempting to initiate a least squares solution ... " << endl;	
		SquareMatrix ATA = A.LeastSquaresHelper(b,"A'*A");
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}
	
	cout << "Case 5:Testing error message #2:" << endl;
	cout << endl;

	{
		double matrix[8] = {1, 1, 1, 1, 1, 2, 3, 4};  
		Matrix A = Matrix::Test(matrix, 4, 2);
		cout << "The m-by-n non-square matrix of coefficients, A = " << endl;
		cout << A << endl;
		double vector1[5] = {6, 5, 7, 10, 11};  
		Matrix b = Matrix::Test(vector1, 5, 1);
		cout << "The p-by-1 column vector, b  = " << endl;
		cout << b << endl;
		cout << "Attempting to initiate a least squares solution ... " << endl;	
		cout << endl;	
		SquareMatrix ATA = A.LeastSquaresHelper(b,"A'*A");
		cout << "Press return to continue ..." << flush;
		system("read");
		cout << endl;
	}	


	return;
}

int main () {

	for (;;) {		
		cout << endl;		
		cout << "Choose to test one of the following:" << endl;
		cout << endl;
		cout << "  Enter \'A\' for the Ones, Eye and other Constructors Testing" << endl;
		cout << "  Enter \'B\' for the Toeplitz function Testing" << endl;
		cout << "  Enter \'C\' for the Transpose function Testing" << endl;
		cout << "  Enter \'D\' for the Triangular Extraction Testing" << endl;
		cout << "  Enter \'E\' for the LU Decomposition Testing - Algorithm 2.1" << endl;
		cout << "  Enter \'F\' for the LU Decomposition Testing - Algorithm 2.2" << endl;
		cout << "  Enter \'G\' for the Forward Substitution Testing - Algorithm 3.1" << endl;
		cout << "  Enter \'H\' for the Backward Substitution Testing - Algorithm 3.2" << endl;
		cout << "  Enter \'I\' for Solving Systems of Linear Equations in which the number of" << endl;
		cout << "            equations is equal to the number of unknowns" << endl;
		cout << "  Enter \'J\' for Solving Systems of Linear Equations in which the number of" << endl;  
		cout << "            equations is greater than the number of unknowns" << endl;
		cout << endl;		
		cout << ">> ";
		char choice;
		cin >> choice;
		cout << endl;
		switch (choice) {
			case 'A':
			case 'a':	ConstructorsTesting();
						break;
			case 'B':
			case 'b':	ToeplitzTesting();
						break;
			case 'C':
			case 'c':	TransposeTesting();
						break;
			case 'D':
			case 'd':	TriangularExtractionTesting();
						break;
			case 'E':
			case 'e':	LUDecompositionTestingOne();
						break;
			case 'F':
			case 'f':	LUDecompositionTestingTwo();
						break;
			case 'G':
			case 'g':	ForwardSubstitution();
						break;
			case 'H':
			case 'h':	BackwardSubstitution();
						break;
			case 'I':
			case 'i':	UniqueSolution();
						break;
			case 'J':
			case 'j':	LeastSquareSolution();
						break;
		}
		cout << "Enter \'0\' to exit or \'1\' to choose another test" << endl;
		cout << endl;
		cout << ">> ";
		cin >> choice;
		if (choice == '0') {
			cout << endl;
			cout << "Goodbye!" << endl;
			cout << endl;
			return 0;
		}
	}

}
