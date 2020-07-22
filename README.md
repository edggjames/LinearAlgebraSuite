# LinearAlgebraSuite

This folder contains code that I have written to implement an object oriented library for linear algebra in C++. It is divided into three
parts:
- Part 1 implements the Matrix class that provides the support for arbitrarily sized 2-D matrices, a.k.a. m-by-n matrices (see Written_Report_1.pdf).
- Part 2 implements the	SquareMatrix class that extends the Matrix class to provide specialised support for square matrices, a.k.a. n-by-n matrices. This is used to  
implement LU decomposition of a square matrix, using both Gaussian elimination without pivoting, and the more numerically stable Gaussian elimination with partial pivoting.
- Part 3 applies the implementation to several important engineering problems. For example,  
in a system of linear equations where the number of equations is the same as the number  
of unknowns (i.e. the solution is uniquely determined), forward substitution and  
backward substitution algorithms are implemented to reach a solution. In a system of  
equations where the number of equations is greater than the number of unknowns (i.e.  
no exact solution exists), a least squares solver is also implemented.    

Written_Report_1.pdf refers to Part 1.
Written_Report_2.pdf refers to Parts 2 and 3.

NB - I did not author MatrixTest.cpp, and SquareMatrixTest.cpp was adapted from this file.
