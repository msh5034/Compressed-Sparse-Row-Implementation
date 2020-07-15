/*
Project 1
DSA 5005
Matt Huson
7/3/19
This program takes input from a file, and stores it as a "Compressed Sparse
Row" format matrix, using three arrays, followed by displaying the matrix info
and the matrix itself. It then adds two matrices together and displays the result.
*/

#include <iostream>
using namespace std;

class CSR
{
protected:
	int noRows; //Number of rows of the original matrix
	int noCols; //Number of columns of the original matrix
	int noNonSparseValues;//Number of NZVs of the original matrix
	int* valueArray; //array of nonsparse values
	int* IA; //the cumulative number of nonsparse values prior to index row
	int* JA; //column numbers of each nonsparse value

public:
	CSR(); //default constructor
	CSR(int n, int m, int numNZV); //overloaded constructor
	void setArrays(int* IA, int* JA, int* VA);
	CSR* Add(CSR& M); //Matrix Add
	void display_matrix(); //display in matrix form
	void display_valueArray(); //display values
	void display_JA(); //display column numbers of values
	void display_IA(); //display cumulative row value count
	int getRows();
	int getCols();
	int getnumNZV();
	~CSR();//destructor
	friend CSR* operator+(const CSR& A, const CSR& B);//overloaded + operator
	friend ostream& operator<<(ostream& out, CSR& X);//overloaded << operator

};

CSR::~CSR() {
	delete valueArray;
	delete IA;
	delete JA;
}//End destructor

CSR::CSR()
{
	noRows = 0;
	noCols = 0;
	noNonSparseValues = 0;
	//I felt like it was unnecessary/unwise to initialize the arrays here
}//End default constructor


CSR::CSR(int n, int m, int numNZV)
{
	//initialized the arrays using the setArrays method rather than the constructor
	noRows = n;
	noCols = m;
	noNonSparseValues = numNZV;
}//End constructor

void CSR::setArrays(int* IA, int* JA, int* VA) {
	//each array is determined either by reading in the matrix or performing the addition function
	//and then the arrays are initialized within the object using this method
	this->valueArray = VA;
	this->IA = IA;
	this->JA = JA;
}//End setArrays


CSR* CSR::Add(CSR& M) {
	int valsEncounteredA = 0;//values encountered so far during matrix A iteration
	int valsEncounteredB = 0;//values encountered so far during matrix B iteration
	int rowValsA = 0;//number of values on the current row for A
	int rowValsB = 0;//number of values on the current row for B
	int lengthJA = (this->noNonSparseValues + M.noNonSparseValues);//max length of matrixC->JA


	int* newIA = new int[this->noRows];//create IA for matrixC, length will always = noRows
	newIA[0] = 0;//will always be 0 due to nature of our IA
	/* for my first step here, I determined the length of matrixC->JA, while also filling in the values
	of matrixC->IA. I started with the max possible values for the length of JA and the greatest possible
	values for each index of IA, and then as I found common values between the matrices, I decremented
	both the length of JA and the value at the index of IA. I later realized that I had self-imposed some
	unnecessary constraints on myself, and that this loop could be replaced by iterating through a
	temporary 2D array.*/
	for (int i = 0; i < this->noRows; ++i) {
		int indexA;
		int indexB;
		int maxA;
		int maxB;
		if (i < this->noRows - 1) {//since IA doesn't store a value for the final row, the for loop needed
									//be split
			/*max value for this index of matrixC->IA. The max difference between IA indices is the sum of
			//the differences of the two constituent IA indices, so the max value it can have is their
			differences + the previous row.*/
			newIA[i + 1] = (this->IA[i + 1] - this->IA[i]) + (M.IA[i + 1] - M.IA[i]) + newIA[i];
			indexA = this->IA[i];//number of cumulative values after previous row 
			indexB = M.IA[i];//number of cumulative values after previous row
			maxA = this->IA[i + 1];//number of cumulative values after this row
			maxB = M.IA[i + 1];//number of cumulative values after this row

			//This will iterate through matrixA->JA for the values in the current row of the matrix
			for (int j = indexA; j < maxA; j++)
			{
				//This will iterate through matrixB->JA for the values in the current row of the matrix
				for (int k = indexB; k < maxB; k++)
				{
					//if a common column value is found, the length of matrixC->JA is reduced by 1, and
					//so is the value at the next index of IA
					if (this->JA[j] == M.JA[k])
					{
						lengthJA--;
						newIA[i + 1]--;
					}//End if
				}//End for
			}//End for
		}//End if
		else {//the iteration for the final row of the matrices
			indexA = this->IA[i];//number of cumulative values before final row
			indexB = M.IA[i];//number of cumulative values before final row
			maxA = this->noNonSparseValues;//number of cumulative values after final row
			maxB = M.noNonSparseValues;//number of cumulative values after final row

			//This will iterate through matrixA->JA for the values in the final row of the matrix
			for (int j = indexA; j < maxA; j++)
			{
				//This will iterate through matrixB->JA for the values in the final row of the matrix
				for (int k = indexB; k < maxB; k++)
				{	//if a common column value is found, the length of matrixC->JA is reduced by 1, but
					//there is no IA value for this row, so there's nothing to change
					if (this->JA[j] == M.JA[k])
					{
						lengthJA--;
					}//End if
				}//End for
			}//End for
		}//End else
	}//End for
	//JA and valueArray are initialized for matrixC, with (lengthJA) indices:
	int* newJA = new int[lengthJA];
	int* newVA = new int[lengthJA];

	//reset all counts to zero and add a counter for values in matrixC:
	valsEncounteredA = 0;
	valsEncounteredB = 0;
	int totalValsEncountered = 0;

	bool encounteredA = false;//flag for encountering a value in A
	for (int row = 0; row < noRows; ++row) {//loop through rows
		for (int col = 0; col < noCols; ++col) {//loop through columns
			if ((row == noRows - 1 || valsEncounteredA < this->IA[row + 1]) && col == this->JA[valsEncounteredA]) {
				//Encounter value for A
				encounteredA = true;//Set value A flag
				newJA[totalValsEncountered] = col;
				newVA[totalValsEncountered++] = this->valueArray[valsEncounteredA++];
			}//End if
			if ((row == noRows - 1 || valsEncounteredB < M.IA[row + 1]) && col == M.JA[valsEncounteredB]) {
				//Encounter value for B
				if (!encounteredA) {//Test value A flag
					//Encounter value for B, but not A
					newJA[totalValsEncountered] = col;
					newVA[totalValsEncountered++] = M.valueArray[valsEncounteredB++];
				}//End if
				else {
					//Encounter value for B and A
					newVA[totalValsEncountered - 1] += M.valueArray[valsEncounteredB++];
					//since totalVals has been incremented, we need to go back to the previous index in order to
					//perform our addition function
				}//End else
			}//End if
			encounteredA = false;//Reset flag
		}//End for
	}//End for

	CSR* returnCSR = new CSR(this->noRows, this->noCols, lengthJA);//create a CSR to return and initialize the first three fields
	returnCSR->setArrays(newIA, newJA, newVA);//intitialize arrays for the returnCSR
	return returnCSR;

}//End Add

CSR* operator+ (CSR& A, CSR& B)
{
	return A.Add(B);
}//End operator+


void CSR::display_matrix()
{
	int nzvVals = 0;//counter for number of NZV values displayed

	for (int i = 0; i < noRows; i++) {//Iterate through IA/rows
		for (int j = 0; j < noCols; j++) {//Iterate through columns
			if (i < noRows - 1 && IA[i + 1] > nzvVals && JA[nzvVals] == j) {//test if there is a JA value corresponding to 
																			//a column, given an IA value
				cout << valueArray[nzvVals++] << " ";
			}//End if
			/*The below is same test as above, but for the last row of matrix. i == noRows - 1 because that's the index
			of the last row and sub in noNonSparseValues for IA[i + 1], since that index doesn't exist and now we're
			just comparing to the total number of non sparse values*/
			else if (i == noRows - 1 && noNonSparseValues > nzvVals && JA[nzvVals] == j) {
				cout << valueArray[nzvVals++] << " ";
			}//End else if
			else {//output a 0 if IA doesn't increase from one row to the next, or if there is not a JA value corresponding to a column
				cout << "0 ";
			}//End else
		}//End for
		cout << endl;
	}//End for

}//End display_matrix

void CSR::display_valueArray()
{
	for (int i = 0; i < noNonSparseValues; i++)
	{
		cout << valueArray[i] << " ";//output every index of valueArray
	}//End for
	cout << endl;
}//End display_valueArray

void CSR::display_JA()
{
	for (int i = 0; i < noNonSparseValues; i++)
	{
		cout << JA[i] << " ";//output every index of JA
	}//End for
	cout << endl;
}//End display_JA

void CSR::display_IA()
{
	for (int i = 0; i < noRows; i++)
	{
		cout << IA[i] << " ";//output every index of IA
	}//End for
	cout << endl;
}//End display_IA
int CSR::getRows()
{
	return noRows;
}//End getRows
int CSR::getCols()
{
	return noCols;
}//End getCols
int CSR::getnumNZV()
{
	return noNonSparseValues;
}//End getnumNZV

ostream& operator << (ostream& out, CSR& X)
{
	//this is written the same as the display_matrix method, with the only differences being
	//the way argument matrix is handled and sending it to "out" instead of "cout"
	int nzvVals = 0;

	for (int i = 0; i < X.noRows; i++) {
		for (int j = 0; j < X.noCols; j++) {
			if (i < X.noRows - 1 && X.IA[i + 1] > nzvVals && X.JA[nzvVals] == j) {
				out << X.valueArray[nzvVals++] << " ";
			}//End if
			else if (i == X.noRows - 1 && X.noNonSparseValues > nzvVals && X.JA[nzvVals] == j) {
				out << X.valueArray[nzvVals++] << " ";
			}//End else if
			else {
				out << "0 ";
			}//End else
		}//End for
		out << endl;
	}//End for
	return out;
}//End operator <<

int main()
{
	int n, m, numNZV;//Declare CSR constructor inputs
	int testValue, count = 0;//Declare variables for use in reading in matrices
	cin >> n >> m >> numNZV;//Read in constructor inputs
	CSR* matrixA = new CSR(n, m, numNZV);//Declare matrixA and initialize with constructor values
	cout << "MATRIX A ---- Rows, Cols and number of non common values: " << matrixA->getRows() << ", "
		<< matrixA->getCols() << ", " << matrixA->getnumNZV() << endl;

	int* IA = new int[n];//Initialize arrays to correspond to arrays in CSR. These arrays will
	int* JA = new int[numNZV];//receive the read in values. The lengths of the arrays are determined
	int* VA = new int[numNZV];//by the size and contents of the matrix.

	IA[0] = 0;//IA[0] is always 0
	for (int i = 0; i < n; i++) { //Iterate through the rows
		for (int j = 0; j < m; j++) { //Iterate through the columns for a given row
			cin >> testValue;//Read in each value of the matrix
			if (testValue != 0) {//Test to make sure it's a non sparse value
				//If non-sparse:
				VA[count] = testValue;//place the value at the [count] index in the VA array
				JA[count] = j;//place the column number at the [count] index in the JA array
				count++;//increment count, basically moving to the next index in both arrays
			}//End if
		}//End for
		if (i < n - 1) {
			IA[i + 1] = count;
		}//End if
	}//End if
	matrixA->setArrays(IA, JA, VA);
	cout << "The valuesArray for matrix A are : ";
	matrixA->display_valueArray();
	cout << "The JA for matrix A are : ";
	matrixA->display_JA();
	cout << "The IA for matrix A are : ";
	matrixA->display_IA();
	cout << "The matrix A is : " << endl;
	cout << *matrixA;


	count = 0;//Reset count, it is unnecessary to reset testValue as it is always overwritten

	cin >> n >> m >> numNZV;//Read in matrix identifiers
	CSR* matrixB = new CSR(n, m, numNZV);//Dynamically create CSR matrixB
	cout << endl;
	cout << "MATRIX B ---- Rows, Cols and number of non common values: " << matrixB->getRows() << ", "
		<< matrixB->getCols() << ", " << matrixB->getnumNZV() << endl;

	IA = new int[n];//Update arrays using new values read in from file
	JA = new int[numNZV];
	VA = new int[numNZV];

	//This process is the same as for matrixA
	IA[0] = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			cin >> testValue;
			if (testValue != 0) {
				VA[count] = testValue;
				JA[count] = j;
				count++;
			}//End if
		}//End for
		if (i < n - 1) {
			IA[i + 1] = count;
		}//End if
	}//End for

	matrixB->setArrays(IA, JA, VA);//Initialize arrays in matrixB to appropriate values

	cout << "The valuesArray for matrix B are : ";
	matrixB->display_valueArray();
	cout << "The JA for matrix B are : ";
	matrixB->display_JA();
	cout << "The IA for matrix B are : ";
	matrixB->display_IA();
	cout << "The matrix B is : " << endl;
	cout << *matrixB;


	CSR* matrixC = new CSR(n, m, numNZV);//Dynamically create CSR matrixC
	matrixC = (*matrixA) + (*matrixB);//Call overloaded + operator

	cout << endl;
	cout << "MATRIX C = A + B" << endl;
	cout << "The valuesArray for matrix C are : ";
	matrixC->display_valueArray();
	cout << "The JA for matrix C are : ";
	matrixC->display_JA();
	cout << "The IA for matrix C are : ";
	matrixC->display_IA();
	cout << "The matrix C is : " << endl;
	cout << *matrixC;//Call overloaded << operator


	delete(matrixC);//Call destructors
	delete(matrixB);
	delete(matrixA);
	return 0;
}//End main
