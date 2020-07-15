Project 1
============

Matthew Huson
--------------
dsa5005, Summer 2019
-------------------

This program reads in a sparse matrix (where sparse values are always 0) and stores it in 3 one-dimensional arrays.

### Assumptions/Bugs ###

* only works if sparse values are 0
* requires number of columns, rows and number of nonzero values to be read in with matrix

### Classes ###

#### CSR ####

* Class Functions (all public)
  * CSR() - default constructor, sets number of rows, columns and values to 0
  * CSR(n,m,numNZV) - overloaded constructor, sets number of rows, columns and values to n, m and numNZV respectively
  * setArrays(IA,JA,VA) - helper function to populate the array attributes
  * Add(M) - function that adds this matrix to matrix M
  * operator+ - overloaded + operator that performs the Add function
  * display_matrix() - display CSR in sparse format
  * display_valueArray() - display valueArray attribute
  * display_IA() - display cumulative row (IA) attribute
  * display_JA() - display column (JA) attribute
  * getRows(), getColumns() and getnumNZV() - get numbers of each attribute
  * ~CSR() - destructor, deletes dynamic array attributes

### Main Function ###

* reads in number of rows, columns and nonzero values
* initializes CSR instance using the read-in values, and displays the values
* iterates through the matrix using the number of columns, rows and nonzero values as a guide
* populates a row array, a column array and value array
* sets CSR attributes to the populated arrays
* reads in a second matrix using the above methodology
* creates a third matrix by adding the two existing matrices together
* displays the array attributes for the third matrix, along with the matrix in sparse format
* calls destructors for all three matrices

