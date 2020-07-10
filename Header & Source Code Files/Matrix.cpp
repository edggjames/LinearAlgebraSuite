/*
 * Matrix.cpp
 *
 *  Created on: 21.12.17
 *      Author: Edward James - Medical Imaging CDT MRes, 2017-2018. 
 */

// Source file for class Matrix, to define implementations for functions declared in Matrix.h

#include "Matrix.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>  //for abs function
#include <iomanip> //for printing precision

using namespace std;

// Part 1.

//getter - to specify the position (index) along 'data' of a matrix entry in the row and column specified by 'rowIdx' and 'columnIdx', respectively
//assuming that the matrix indexing starts at (0,0) and the array indexing starts at (0).
int Matrix::GetIndex (const int& rowIdx, const int& columnIdx) const {
	int index = rowIdx + columnIdx*noOfRows; 
	return index;
}

// 1. ConstructorsTesting

//default constructor
Matrix::Matrix () {
	noOfRows = 0;
	noOfCols = 0;
	data = 0;
}

//standard custom constructor
Matrix::Matrix (const int& noOfRows_rhs, const int& noOfCols_rhs) {
	// first check that both inputs are greater than zero
	if (noOfRows_rhs <= 0 or noOfCols_rhs <= 0) {
		cerr << endl;		
		cerr << "Error:  In standard custom constructor ..." << endl;
		cerr << "	 Cannot create a matrix with a non-positive dimension." << endl;
		cerr << "        Exiting program ... " << endl;
		cerr << endl;
		exit(1);
	}
	//assign dimension values
	noOfRows = noOfRows_rhs;
	noOfCols = noOfCols_rhs;	
	//allocate new memory
	int product = noOfRows*noOfCols;
	data = new double [product];
	//then assign '0' to each value (this is done linearly here, rather than column-wise, as all values are the same)
	for (int i = 0; i < product; ++i) {
		data[i] = 0;
	}
}

// default copy constructor
Matrix::Matrix (const Matrix& rhs) {
	
	//assign dimensions	
	this->noOfRows = rhs.noOfRows;
	this->noOfCols = rhs.noOfCols;

	// allocate the memory of appropriate size
	int product = noOfRows*noOfCols;	
	this->data = new double[product];

	// ensure that a 'deep', rather than a 'shallow', copy is done
	for (int i = 0; i < product; ++i) {
		this->data[i] = rhs.data[i];	
	}		
}

// static Zeros constructor
Matrix Matrix::Zeros(const int& noOfRows_rhs, const int& noOfCols_rhs) {
	// There is no vector object attached to a static member function
	
	// Therefore need to create Matrix object of zeros using standard custom constructor
	Matrix Zeros (noOfRows_rhs, noOfCols_rhs);	

	// return the object
	return Zeros;
}

// static Ones constructor
Matrix Matrix::Ones(const int& noOfRows_rhs, const int& noOfCols_rhs) {
	
	//Create Matrix object of zeros using standard custom constructor
	Matrix Ones (noOfRows_rhs, noOfCols_rhs);	

	//Then assign '1' to each value (this is done linearly here, rather than column-wise, as all values are the same)
	int product = noOfRows_rhs*noOfCols_rhs;	
	for (int i = 0; i < product; ++i) {
		Ones.data[i] = 1;
	}
	
	//Return the object
	return Ones;
}


//destructor 
Matrix::~Matrix () {
	//dynamically delete all current instances of member variable data
	delete[] data;
}

//printer (remembering that the matrix has been stored columnwise in data)
ostream& operator<< (ostream& out, const Matrix& rhs)  {			
	//loop through the row indices
	for (int i = 0; i < rhs.noOfRows; ++i) { 
	out << "\t"; // leave a tab (8 spaces) at the beginning of the row
		//loop through the column indices
		for (int j = 0; j < rhs.noOfCols; ++j) {			
			int index = i + j*rhs.noOfRows;		// GetIndex (i,j) cannot be used here as not declared in this scope	
			double element = rhs.data[index];
			double tolerance = 1e-14;
			if (abs(element) < tolerance) { //set to zero if magnitude of element is less than tolerance;
				element = 0;
			}			
			if (element >= 0) {		
				out << " " << element; // leave an extra space if element is positive, to ensure alignment of matrix entries in console out
			}
			else {
				out << element;
			}
			out << "         "; // leave 9 spaces between columns (tab + 1 doesn't work ...)
		}
		out << endl; // go to next line at end of row
	}
	
	return out;	
}

// 2. AssignmentTesting

// default assignment operator
Matrix& Matrix::operator= (const Matrix& rhs) {
	// check the dimensions of the rhs input
	if (this->noOfRows != rhs.noOfRows and this->noOfCols!= rhs.noOfCols) {
		cerr << "Error:  Assignment operation failed ..." << endl;
		cerr << "        rhs input has incompatible number of rows and columns." << endl;
		cerr << endl;
	}	
	else if (this->noOfRows != rhs.noOfRows) {
		cerr << "Error:  Assignment operation failed ..." << endl;
		cerr << "        rhs input has incompatible number of rows." << endl;
		cerr << endl;	
	}
	else if (this->noOfCols!= rhs.noOfCols) {
		cerr << "Error:  Assignment operation failed..." << endl;
		cerr << "        rhs input has incompatible number of columns." << endl;
		cerr << endl;	
	}
	//guard against copy into self
	else if (this == &rhs) { 
		cerr << "Error:  Assignment operation failed..." << endl;
		cerr << "        copy into self attempted." << endl;
		cerr << endl;
	}
	else {
		// assign the data of rhs to lhs
		int product = noOfRows*noOfCols;
		for (int i = 0; i < product; ++i) {		
			this->data[i] = rhs.data[i];	
		}
	}
	// return the reference to self
	return *this;	
}

// 3. ToeplitzTesting & ToeplitzTestingHelper

//print function
Matrix Matrix::Print (const double* const data_rhs, const int& noOfRows_rhs, const int& noOfCols_rhs) {	
	
	//Create Matrix object of zeros using standard custom constructor
	Matrix Print (noOfRows_rhs, noOfCols_rhs);

	//copy across the content of data_rhs
	int product = Print.noOfRows*Print.noOfCols;
	for (int i = 0; i < product; ++i) {		
			Print.data[i] = data_rhs[i];	
		}

	//loop through the row indices of print
	for (int i = 0; i < Print.noOfRows; ++i) { 
	cout << "\t"; // leave a tab (8 spaces) at the beginning of the row	
		//loop through the column indices of print
		for (int j = 0; j < Print.noOfCols; ++j) {			
			int index = i + j*Print.noOfRows;   	// GetIndex(i,j) cannot be user here without object	
			double element = data_rhs[index];
			if (element >= 0) {		
				cout << " " << element; // leave an extra space if element is positive, to ensure alignment of matrix entries in console output
			}
			else {
				cout << element;
			}
			cout << "         "; // leave 9 spaces between columns (tab + 1 doesn't work ...)
		}
		cout << endl; // go to next line at end of row
	}
	
	return Print;
}

//generate Toeplitz function
Matrix Matrix::Toeplitz (const double* const column, const int& noOfRows_rhs, const double* const row, const int& noOfCols_rhs) {
	
	//Create Matrix object of zeros using standard custom constructor
	Matrix Toeplitz (noOfRows_rhs, noOfCols_rhs);

	//copy across the content of data_rhs
	int product = Toeplitz.noOfRows*Toeplitz.noOfCols;
	
	//Display a warning message if first element of column and first element of row are not equal
	if (column[0] != row[0]) {
		cout << "Warning: In Toeplitz constructor, first element of input column does not match" << endl;
		cout << "         first element of input row! Column wins diagonal conflict." << endl;
		cout << endl;
	}

	//To construct Toeplitz matrix based on columns and row, compressed into a 1D array columnwise
	//loop through the rows indices
	for (int i = 0 ; i <Toeplitz.noOfRows; i++) {
		//loop through the column indices
		for (int j = 0 ; j <Toeplitz.noOfCols; j++) {	
			int index = (i + j*Toeplitz.noOfRows);	// GetIndex(i,j) cannot be user here without object
			//for diagonal of matrix or below the diagonal	
			if (i>=j){			
				Toeplitz.data[index]=column[i-j];
			}
			//above this diagonal line
			else if (j>i){
				Toeplitz.data[index]=row[j-i];
			}
		}
	}

	// return the object
	return Toeplitz;	
}

// 4. TransposeTesting

//static transpose function
Matrix Matrix::Transpose(const Matrix& rhs) {
	//create an instance of Matrix object called Transpose
	Matrix Transpose(rhs); 	//calls default copy constructor to replicate rhs	
	
	//invoke the non-static transpose function (see below)	
	Transpose.Transpose();
	
	//return the object
	return Transpose;

}

//non-static transpose function
void Matrix::Transpose () {
	//swap the number of columns and rows
	int tmp = 0;	
	tmp = noOfRows;
	noOfRows = noOfCols;
	noOfCols = tmp;

	// allocate the memory for temporary holder for data
	int product = noOfRows*noOfCols;
	double *tmp_data = new double[product];
	
	//re-allocate members of transposed data as a tranpose of original data
	int k = 0; //use this as a counter to index original data
	//loop through the row indices of tranposed data
	for (int i = 0; i < noOfRows; ++i) {
		//looping through the column indices of transposed data
		for (int j = 0; j < noOfCols; ++j) {
			int index = GetIndex(i,j); //to index tranposed data
			tmp_data[index] = data[k]; //appropriate transformation to switch columns and rows (hold temporarily in tmp_data)	
			k += 1;	//increment the counter			
			}
		}
	
	// then copy across the content from tmp_data to data in a linear fashion
	for (int i = 0; i < product; ++i) {		
			data[i] = tmp_data[i];	
	}
	
	//delete dynamically allocated memory for tmp_data
	delete [] tmp_data;

	//empty return as void function
	return;
}

// 5. MultiplicationTesting

//non-member multiplication function - creates new matrix
Matrix operator* (const Matrix& lhs, const Matrix& rhs) {	
	
	Matrix result(lhs); //calls default copy constructor
	result *= rhs;	    //calls self-modifying member multiplication function (see below)
	return result;	    //returns object
	
}

//member multiplication function
Matrix& Matrix::operator*= (const Matrix& rhs) {	
	
	//check that number of cols of LHS is equal to number of rows of RHS
	if (this->noOfCols != rhs.noOfRows) {
		cerr << endl;
		cerr << "Error:  In member multiplication function ... " << endl;
		cerr << "        Inner matrix dimensions of rhs and lhs do not agree." << endl;
		cerr << "        Unmultiplied lhs matrix has been returned." << endl;
		cerr << endl;
		Matrix error;
		return *this;
	}
	// The returned matrix will have the same number of rows as LHS matrix
	// But the returned matrix will have the same number of columns as the RHS matrix
	// Therefore amend number of columns of returned matrix if need be
	if (this->noOfCols != rhs.noOfCols) {
		int noOfCols_old = this->noOfCols;		
		int product_old = this->noOfRows*this->noOfCols;		
		this->noOfCols = rhs.noOfCols;		
		int product_new = this->noOfRows*this->noOfCols;		
		if (noOfCols_old < this->noOfCols) {
			//extend data by allocating extra columns and padding with zeros:		
			//dynamically allocate memory for array of zeros of appropriate length
			double *tmp = new double[product_new];
			//initialise these elements all to zero
			for (int i = 0; i < product_new; ++i) {
				tmp[i] = 0;
			}
			//then copy across the content from data to temp
			for (int i = 0; i < product_old; ++i) {		
				tmp[i] = this->data[i];
			}
			//reassign data 
			this->data = tmp; 
 		}	
	}

	// allocate the memory for temporary holder for data
	int product = this->noOfRows*this->noOfCols;
	double *tmp_data = new double[product];
	 	
	// looping through the row indices of lhs	
	for (int i = 0; i < this->noOfRows; ++i) {
		//looping through the column indices of lhs
		for (int j = 0; j < this->noOfCols; ++j) {					
			double tmp = 0; //to store temporary element value of tmp_data in
        		//looping through the rows of the rhs matrix
			for (int k = 0; k < rhs.noOfRows; ++k) {
       				int index = GetIndex(i,j); //use to index tmp_data to temporarily store result in
				//tmp_data[index] = dot product of (lhs row i) and (rhs col j)				
				int index_LHS = i + k*this->noOfRows; // use to index lhs matrix for calculation				
				int index_RHS = k + j*rhs.noOfRows; //use to index rhs matrix for calculation 				
				tmp += this->data[index_LHS] * rhs.data[index_RHS];  //increment tmp by this result						
				tmp_data[index] = tmp;	
			}		
		}	
	}

	// then copy across the content from tmp_data to data	
	for (int i = 0; i < product; ++i) {		
			this->data[i] = tmp_data[i];	
	}
	
	//delete dynamically allocated memory for tmp_data
	delete [] tmp_data;

	// return the reference to self
	return *this;
}

// 6. RowColumnExchangeTesting

//ExchangeRows function (overloaded with two input parameter options)

//swaps row1 and row2, across all columns
Matrix& Matrix::ExchangeRows (const int& row1, const int& row2) { 
	
	//looping through all columns
	for (int j = 0; j < noOfCols; ++j) {
		int index_row1 = GetIndex(row1,j);    //get index for element in row1
		double tmp = data[index_row1];        //store in tmp variable
		int index_row2 = GetIndex(row2,j);    //get index for element in row2 
		data[index_row1] = data[index_row2];  //swap the elements
		data[index_row2] = tmp;		      //swap the elements
	}

	// return the reference to self
	return *this;	
}

//swaps row1 and row2, but only between col1 and col2
Matrix& Matrix::ExchangeRows (const int& row1, const int& row2, const int& col1, const int& col2) { 

	//looping through columns from col1 to col2
	for (int j = col1; j <= col2; ++j) {
		int index_row1 = GetIndex(row1,j);    //get index for element in row1
		double tmp = data[index_row1];        //store in tmp variable
		int index_row2 = GetIndex(row2,j);    //get index for element in row2 
		data[index_row1] = data[index_row2];  //swap the elements
		data[index_row2] = tmp;		      //swap the elements
	}

	// return the reference to self
	return *this;	
}

//ExchangeColumns function (overloaded with two input parameter options)

//swaps col1 and col2, across all rows
Matrix& Matrix::ExchangeColumns (const int& col1, const int& col2) {

	//looping through all rows
	for (int i = 0; i < noOfRows; ++i) {
		int index_col1 = GetIndex(i,col1);    //get index for element in col1
		double tmp = data[index_col1];        //store in tmp variable
		int index_col2 = GetIndex(i,col2);    //get index for element in col2 
		data[index_col1] = data[index_col2];  //swap the elements
		data[index_col2] = tmp;		      //swap the elements
	}

	// return the reference to self
	return *this;	

}

//swaps col1 and col2, but only between row1 and row2
Matrix& Matrix::ExchangeColumns (const int& col1, const int& col2, const int& row1, const int& row2){

	//looping through rows between row1 and row2
	for (int i = row1; i <= row2; ++i) {
		int index_col1 = GetIndex(i,col1);    //get index for element in col1
		double tmp = data[index_col1];        //store in tmp variable
		int index_col2 = GetIndex(i,col2);    //get index for element in col2 
		data[index_col1] = data[index_col2];  //swap the elements
		data[index_col2] = tmp;		      //swap the elements
	}

	// return the reference to self
	return *this;	

}

// 7. OtherTesting

//Zeros (sets every entry to zero)
void Matrix::Zeros () {
	//Assign 0 to every element in data, in a linear fashion	
	int product = noOfRows*noOfCols;
	for (int i = 0; i < product; ++i) {
		data[i] = 0;
	}

	//empty return as void function
	return;	
}

//Ones (sets every entry to one)
void Matrix::Ones () {
	//Assign 1 to every element in data, in a linear fashion	
	int product = noOfRows*noOfCols;
	for (int i = 0; i < product; ++i) {
		data[i] = 1;
	}

	//empty return as void function
	return;		
}

//GetNoOfRows (find out number of rows)
int Matrix::GetNoOfRows () const {
	return noOfRows;
}	

//GetNoOfColumns (find out number of columns)
int Matrix::GetNoOfColumns () const {	
	return noOfCols;
}

//GetEntry (find out value of particular entry (i,j) in Matrix, starting indexing at (0,0) in upper left corner of matrix)
double Matrix::GetEntry (const int& rowIdx, const int& columnIdx) const {
	int index= GetIndex (rowIdx, columnIdx);
	return data[index];
}

// Part 2.

//To construct user specified matrices for testing functions
Matrix Matrix::Test (const double* data_rhs, const int& noOfRows_rhs, const int& noOfCols_rhs) {

	// create an empty Matrix object
	Matrix Test;
	
	// then set up the dimensions
	Test.noOfRows = noOfRows_rhs;
	Test.noOfCols = noOfCols_rhs;
	
	// allocate the memory
	int product = noOfRows_rhs*noOfCols_rhs;
	Test.data = new double[product];
	
	//then copy across the content from data_rhs to Test.data
	for (int i = 0; i < product; ++i) {
		Test.data[i] = data_rhs[i];
	}
	
	// return the object
	return Test;
}

//Print function to print matrices in default MATLAB style to fixed 4.d.p precision
Matrix Matrix::PrintMATLAB (const Matrix& rhs) {
	//Need to create an empty SquareMatrix object 
	Matrix print;

	// then set up the dimensions
	print.noOfRows = rhs.noOfRows;
	print.noOfCols = rhs.noOfCols;

	// allocate the memory
	int product = print.noOfRows*print.noOfCols;
	print.data = new double[product];

	//copy across the content of data_rhs
	for (int i = 0; i < product; ++i) {		
			print.data[i] = rhs.data[i];	
		}

	//set to display numbers to fixed 4.d.p.
	cout << fixed << setprecision(4);

	//loop through the row indices of print
	for (int i = 0; i < print.noOfRows; ++i) { 
	cout << "\t"; // leave a tab (8 spaces) at the beginning of the row	
		//loop through the column indices of print
		for (int j = 0; j < print.noOfCols; ++j) {			
			int index = i + j*print.noOfRows;   	// GetIndex(i,j) cannot be user here without object	
			double element = print.data[index];
			double tolerance = 1e-14;		
			if (abs(element) < tolerance) {
				element = 0;  //set to zero if magnitude of element is less than tolerance;
			}  
			if (element == 0) {
				cout << "      0";  //leave 6 spaces to ensure alignment
			}			
			else if (element >= 0) {		
				cout << " " << element; // leave an extra space if element is positive, to ensure alignment of matrix entries in console output			
			}
			else {
				cout << element;
			}
			cout << "         "; // leave 9 spaces between columns (tab + 1 doesn't work ...)
		}
		cout << endl; // go to next line at end of row
	}
	cout << endl;
		
	cout.unsetf(ios::floatfield);       // unset floatfield from 'fixed' to 'not set'
  	cout << setprecision(6);  	    // set precision back to default value of 6.d.p

	return print;
}

// Part 3. 

// 8. Forward Substitution - Algorithm 3.1 - implemented here so can access and manipulate column vector b

Matrix Matrix::ForwardSub (const Matrix& b) const {
	
	//construct a m-by-1 column vector to hold return result y in
	int dim = this->noOfCols;	
	Matrix y = Matrix(dim,1);

	//Expanding from the top row of L:
	//1. Calculate first entry of y
	y.data[0] = b.data[0]/this->data[0];
	
	//2. Calculate subsequent entries of y	
	for (int k= 1; k <dim; ++k) {
		double tmp = 0;		
		for (int i= 0; i <k; ++i) {
			int index_ki = GetIndex(k,i);
			tmp += (this->data[index_ki]*y.data[i]);	
		}		
		int index_kk = GetIndex(k,k);
		y.data[k] =  (b.data[k] - tmp)/this->data[index_kk]; 
	
	}

	return y;

}

// 9. Backward Substitution - Algorithm 3.2 - implemented here so can access and manipulate column vector y

Matrix Matrix::BackwardSub (const Matrix& y) const {

	//construct a m-by-1 column vector to hold return result x in
	int dim = this->noOfCols;
	Matrix x = Matrix(dim,1);

	//Expanding from the bottom row of U:
	//1. Calculate first entry of x
	x.data[dim-1] = y.data[dim-1]/this->data[dim*dim-1];
	
	//2. Calculate subsequent entries of x
	for (int k= dim-2; k >=0; --k) {
		//cout << k << endl;
		double tmp = 0;		
		for (int i= k+1; i <dim; ++i) {
			int index_ki = GetIndex(k,i);
			tmp += (this->data[index_ki]*x.data[i]);	
		}		
		int index_kk = GetIndex(k,k);
		x.data[k] =  (y.data[k] - tmp)/this->data[index_kk]; 
	}

	return x;

}

// 10. Solving a system of linear equations in which the number of equations is equal to the number of unknowns, i.e. a unique determined solution exists. 

Matrix Matrix::SolveUnique (const Matrix& L, const Matrix& U, const Matrix& P, Matrix& b) const {
	
	//multiply P by b
	b = P*b;  		 //NB assignment operator throws an error here if number of columns in b and P do not match

	//solve for y in Ly=b using Algorithm 3.1
	Matrix y = L.ForwardSub(b); 
	
	//solve for x in Ux=y using Algorithm 3.2
	Matrix x = U.BackwardSub(y);

	return x;
} 

// 11. Solving an overdetermined system of linear equations in which the number of equations is greater than the number of unknowns, i.e. a least squares solution is produced. 

Matrix Matrix::LeastSquaresHelper (const Matrix& b, const string& output) const {

	//check that number of rows in A is bigger than the number of columns in A	
	if (this->noOfRows == this->noOfCols and output == "A'*A") { 
		cout << "Warning: In function LeastSquaresHelper ..." << endl;
		cout << "         Matrix A is square." << endl;
		cout << "         An exact solution, not a least squares solution, will be provided." << endl;
		cout << "         Consider using SolveUnique instead." << endl; 
		cout << endl;
	}
	// check that m > n - error #1
	if (this->noOfRows < this->noOfCols) {
		cerr << endl; 
		cerr << "Error:  In function LeastSquaresHelper ... " << endl;
		cerr << "        In matrix A, the number of rows is less than the number of columns." << endl;
		cerr << endl;
		Matrix error;  
		return error;
	}
	//check that b is of appropriate length	-  error #2
	if (this->noOfRows != b.noOfRows) { 
		cerr << "Error:  In function LeastSquaresHelper ... " << endl;
		cerr << "        The number of rows in b does not equal the number of rows in A." << endl;
		cerr << endl;
		Matrix error;  
		return error;
	}

	//produce ATA
	Matrix A (*this);
	Matrix AT = Matrix::Transpose(A);
	Matrix ATA = AT*A;

	//produce ATb
	Matrix ATb = AT*b; 
	
	//return ATA or ATb dependent upon output string
	if (output == "A'*A") { 			
		cout << "Number of rows in A'*A = " << ATA.GetNoOfRows()  << endl;
		cout << "Number of columns in A'*A = " << ATA.GetNoOfColumns() << endl;
		cout << endl; 
		return ATA;
	}
	if (output == "A'*b") {
		return ATb;
	}
}
