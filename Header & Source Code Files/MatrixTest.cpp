	/*
 * MatrixTest.cpp
 *
 *      Author: Gary Hui Zhang (gary.zhang@ucl.ac.uk)
 */

#include "Matrix.h"
#include <iostream>
#include <cstdlib>  //need this for 'system' functions

using namespace std;

void ConstructorsTesting () {
	cout << "Testing the Matrix constructors:" << endl;

	cout << "Case 1: Creating a 2x4 matrix of zeros with the standard constructor:" << endl;
	{
		Matrix matrix(2, 4);
		cout << matrix << endl;		
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: Creating a 2x4 matrix of zeros with the static Zeros:" << endl;
	{
		Matrix matrix = Matrix::Zeros(2, 4);
		cout << matrix << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: Creating a 2x4 matrix of ones with the static Ones:" << endl;
	{
		Matrix matrix = Matrix::Ones(2, 4);
		cout << matrix << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 4: Copying a 2x4 matrix of ones with the copy constructor:" << endl;
	{
		Matrix matrix = Matrix::Ones(2, 4);
		cout << "The input matrix = " << endl;
		cout << matrix << endl;
		Matrix copy = matrix; 
		cout << "The copy = " << endl;
		cout << copy << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}
}

void AssignmentTesting () {
	cout << "Testing the assignment operator:" << endl;

	cout << "Case 1: Using the assignment operator to assign a 2x4 matrix of zeros to a 2x4 matrix of ones:" << endl;
	{
		Matrix zeros = Matrix::Zeros(2, 4);
		cout << "The matrix of zeroes = " << endl;
		cout << zeros << endl;
		Matrix ones = Matrix::Ones(2, 4);		
		cout << "The matrix of ones before assignment = " << endl;
		cout << ones << endl;			
		ones = zeros;
		cout << "The matrix after assignment = " << endl;
		cout << ones << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: Attempting to use the assignment operator to assign a 3x4 matrix of zeros to a 2x4 matrix of ones:" << endl;
	{
		Matrix zeros = Matrix::Zeros(3, 4);
		cout << "The matrix of zeroes = " << endl;
		cout << zeros << endl;
		Matrix ones = Matrix::Ones(2, 4);		
		cout << "The matrix of ones before assignment = " << endl;
		cout << ones << endl;			
		ones = zeros;
		cout << "The matrix after attempted assignment = " << endl;
		cout << ones << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: Attempting to use the assignment operator to assign a 2x3 matrix of zeros to a 2x4 matrix of ones:" << endl;
	{
		Matrix zeros = Matrix::Zeros(2, 3);
		cout << "The matrix of zeroes = " << endl;
		cout << zeros << endl;
		Matrix ones = Matrix::Ones(2, 4);		
		cout << "The matrix of ones before assignment = " << endl;
		cout << ones << endl;			
		ones = zeros;
		cout << "The matrix after attempted assignment = " << endl;
		cout << ones << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}	

	cout << "Case 4: Attempting to use the assignment operator to assign a 3x5 matrix of zeros to a 2x4 matrix of ones:" << endl;
	{
		Matrix zeros = Matrix::Zeros(3, 5);
		cout << "The matrix of zeroes = " << endl;
		cout << zeros << endl;
		Matrix ones = Matrix::Ones(2, 4);		
		cout << "The matrix of ones before assignment = " << endl;
		cout << ones << endl;			
		ones = zeros;
		cout << "The matrix after attempted assignment = " << endl;
		cout << ones << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

}

void ToeplitzTestingHelper (const double* const column, const int& noOfRows, const double* const row, const int& noOfColumns, const double* const expected) {
	cout << "The 1st column = " << endl;
	Matrix::Print(column, noOfRows, 1);
	cout << endl;
	cout << "The 1st row = " << endl;
	Matrix::Print(row, 1, noOfColumns);
	cout << endl;
	cout << "The matrix created by the toeplitz function in MATLAB = " << endl;
	Matrix::Print(expected, noOfRows, noOfColumns);
	cout << endl;
	Matrix toeplitz = Matrix::Toeplitz(column, noOfRows, row, noOfColumns);
	cout << "The matrix created by Matrix::Toeplitz = " << endl;
	cout << toeplitz << endl;
}

void ToeplitzTesting () {
	cout << "Testing the static function Matrix::Toeplitz:" << endl;

	cout << "Case 1: When the number of columns is less than the number of rows" << endl;  //changed a typo here
	{
		double column[4] = {2, 1, 0, -1};
		double row[3] = {2, 0, -1};
		// matrix should be stored, in 1-D, column-wise
		//  2     0    -1
		//  1     2     0
		//  0     1     2
		// -1     0     1
		double expected[12] = {2, 1, 0, -1, 0, 2, 1, 0, -1, 0, 2, 1};
		ToeplitzTestingHelper(column, 4, row, 3, expected);
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: When the number of columns is larger than the number of rows" << endl; //changed a typo here
	{
		double column[3] = {2, 0, -1};
		double row[4] = {2, 1, 0, -1};
		// matrix should be stored, in 1-D, column-wise
		//  2     1     0    -1
		//  0     2     1     0
		// -1     0     2     1
		double expected[12] = {2, 0, -1, 1, 2, 0, 0, 1, 2, -1, 0, 1};
		ToeplitzTestingHelper(column, 3, row, 4, expected);
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 3: When the 1st entries of the 1st column and row are different:" << endl;
	{
		double column[4] = {2, 1, 0, -1};
		double row[3] = {1, 0, -1};
		// matrix should be stored, in 1-D, column-wise
		//  2     0    -1
		//  1     2     0
		//  0     1     2
		// -1     0     1
		double expected[12] = {2, 1, 0, -1, 0, 2, 1, 0, -1, 0, 2, 1};
		ToeplitzTestingHelper(column, 4, row, 3, expected);
	}
}

void TransposeTesting () {
	cout << "Testing the Transpose functions:" << endl;

	cout << "Case 1: the static Transpose function" << endl;
	{
		// the same matrix as in ToeplitzTesting
		double column[4] = {2, 1, 0, -1};
		double row[3] = {2, 0, -1};
		Matrix matrix = Matrix::Toeplitz(column, 4, row, 3);
		cout << "The original Matrix = " << endl;
		cout << matrix << endl;
		Matrix transpose = Matrix::Transpose(matrix);
		cout << "The transposed version = " << endl;		
		cout << transpose << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: the non-static Transpose function" << endl;
	{
		// the same matrix as in ToeplitzTesting
		double column[4] = {2, 1, 0, -1};
		double row[3] = {2, 0, -1};
		Matrix matrix = Matrix::Toeplitz(column, 4, row, 3);
		cout << "The original Matrix = " << endl;
		cout << matrix << endl;
		matrix.Transpose();
		cout << "The transposed version = " << endl;
		cout << matrix << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

}

void MultiplicationTesting () {
	cout << "Testing the multiplication functions:" << endl;

	cout << "Case 1: the non-member multiplication function" << endl;
	{
		// the same matrix as in ToeplitzTesting
		double column[4] = {2, 1, 0, -1};
		double row[3] = {2, 0, -1};
		Matrix lhs = Matrix::Toeplitz(column, 4, row, 3);
		cout << "The L.H.S. Matrix = " << endl;
		cout << lhs << endl;
		Matrix rhs = Matrix::Transpose(lhs);
		cout << "The R.H.S. version = " << endl;
		cout << rhs << endl;
		cout << "The expected product = " << endl;
		//  5     2    -2    -3
		//  2     5     2    -1
		// -2     2     5     2
		// -3    -1     2     2
		double expected[16] = {5, 2, -2, -3, 2, 5, 2, -1, -2, 2, 5, 2, -3, -1, 2, 2};
		Matrix::Print(expected, 4, 4);
		cout << "The actual product = " << endl;
		cout << lhs * rhs << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: the member multiplication function" << endl;
	{
		// the same matrix as in ToeplitzTesting
		double column[4] = {2, 1, 0, -1};
		double row[3] = {2, 0, -1};
		Matrix lhs = Matrix::Toeplitz(column, 4, row, 3);
		cout << "The L.H.S. Matrix = " << endl;
		cout << lhs << endl;
		Matrix rhs = Matrix::Transpose(lhs);
		cout << "The R.H.S. version = " << endl;
		cout << rhs << endl;
		cout << "The expected product = " << endl;
		//  5     2    -2    -3
		//  2     5     2    -1
		// -2     2     5     2
		// -3    -1     2     2
		double expected[16] = {5, 2, -2, -3, 2, 5, 2, -1, -2, 2, 5, 2, -3, -1, 2, 2};
		Matrix::Print(expected, 4, 4);
		cout << "The actual product = " << endl;
		lhs *= rhs;
		cout << lhs << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

}

void RowColumnExchangeTesting () {
	cout << "Testing the Row or Column Exchange functions:" << endl;

	cout << "Case 1: the row exchange function" << endl;
	{
		// the same matrix as in ToeplitzTesting
		double column[4] = {2, 1, 0, -1};
		double row[3] = {2, 0, -1};
		Matrix matrix = Matrix::Toeplitz(column, 4, row, 3);
		cout << "The original Matrix = " << endl;
		cout << matrix << endl;
		cout << "Exchange the 2nd and 4th rows gives" << endl;
		cout << matrix.ExchangeRows(1, 3) << endl;
		cout << "Exchange the 1st and 4th rows, but only from the 2nd to the 3rd columns, gives" << endl;
		cout << matrix.ExchangeRows(0, 3, 1, 2) << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Case 2: the column exchange function" << endl;
	{
		// the same matrix as in ToeplitzTesting
		double column[4] = {2, 1, 0, -1};
		double row[3] = {2, 0, -1};
		Matrix matrix = Matrix::Toeplitz(column, 4, row, 3);
		cout << "The original Matrix = " << endl;
		cout << matrix << endl;
		cout << "Exchange the 1st and 2nd columns gives" << endl;
		cout << matrix.ExchangeColumns(0, 1) << endl;
		cout << "Exchange the 2nd and 3rd columns, but only from the 2nd to the 4th rows, gives" << endl;
		cout << matrix.ExchangeColumns(1, 2, 1, 3) << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

}

void OtherTesting () {
	cout << "Testing the miscellaneous functions:" << endl;

	cout << "The function that sets every entry to zero" << endl;
	{
		Matrix ones = Matrix::Ones(2, 4);
		cout << "The input matrix =" << endl;
		cout << ones << endl;
		ones.Zeros();
		cout << "After setting every entry to zero, we have" << endl;
		cout << ones << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "The function that sets every entry to one" << endl;
	{
		Matrix zeros = Matrix::Zeros(2, 4);
		cout << "The input matrix =" << endl;
		cout << zeros << endl;
		zeros.Ones();
		cout << "After setting every entry to one, we have" << endl;
		cout << zeros << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Find out the number of rows and columns" << endl;
	{
		Matrix zeros = Matrix::Zeros(2, 4);
		cout << "The input matrix =" << endl;
		cout << zeros << endl;
		cout << "The number of rows = " << zeros.GetNoOfRows() << " ; ";
		cout << "the number of columns = " << zeros.GetNoOfColumns() << endl;
		cout << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}

	cout << "Find out the value of a particular entry" << endl;
	{
		// the same matrix as in ToeplitzTesting
		double column[4] = {2, 1, 0, -1};
		double row[3] = {2, 0, -1};
		Matrix matrix = Matrix::Toeplitz(column, 4, row, 3);
		cout << "The original Matrix = " << endl;
		cout << matrix << endl;
		cout << "The entry at the 2nd row and the 3rd column is ";
		cout << matrix.GetEntry(1, 2) << endl;
		cout << "The entry at the 3rd row and the 2nd column is ";
		cout << matrix.GetEntry(2, 1) << endl;
		cout << endl;
		cout << "Press any key to continue ..." << flush;
		system("read");
		cout << endl;
	}
}

int main () {

	for (;;) {
		cout << "Choose to test one of the following:" << endl;
		cout << "  Enter \'1\' for the constructors" << endl;
		cout << "  Enter \'2\' for the assignment operator" << endl;
		cout << "  Enter \'3\' for the Toeplitz function" << endl;
		cout << "  Enter \'4\' for the Transpose function" << endl;
		cout << "  Enter \'5\' for the multiplication" << endl;
		cout << "  Enter \'6\' for the Row or Column Exchange functions" << endl;
		cout << "  Enter \'7\' for the other functions" << endl;
		cout << ">> ";
		char choice;
		cin >> choice;
		switch (choice) {
			case '1':	ConstructorsTesting();
						break;
			case '2':	AssignmentTesting();
						break;
			case '3':	ToeplitzTesting();
						break;
			case '4':	TransposeTesting();
						break;
			case '5':	MultiplicationTesting();
						break;
			case '6':	RowColumnExchangeTesting();
						break;
			case '7':	OtherTesting();
						break;
		}
		cout << "Enter \'0\' to exit or \'1\' to choose another test" << endl;
		cout << ">> ";
		cin >> choice;
		if (choice == '0') {
			return 0;
		}
	}

}
