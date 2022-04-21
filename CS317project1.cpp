// CS317project1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stack>
#include <cmath>

using namespace std;

int comparison_counter = 0; 

bool LESS(double a, double b);

void swap(double *a, double *b);

void QuickSort(double A[], int left, int right, int n);
int Partition(double A[], int left, int right);

void MergeSort(double A[], int n);
void Merge(double B[], double C[], double A[], int p, int q);



int main()
{

	ifstream infile;
	infile.open("input.txt");
	if (!infile) { cout << "Error: The input file could not be opened" << endl; }
	else { cout << "The input file was successfully opened!" << endl; }

	ofstream outfile;
	outfile.open("rmk0017_quick.txt"); cout << "An output file has been created!" << endl;

	unsigned int size; 
	infile >> size;
	cout << "There are " << size << " elements in the array" << endl;

	double *list = new double[size];

	for (int i = 0; i < size; i++) {
		infile >> list[i];
	}

	cout << endl << endl;
	
	cout << "**************************\n";
	cout << "Implementing QuickSort now...\n";
	int left = 0;
	int right = size - 1;
	QuickSort(list, left, right, size);
	cout << "**************************\n";
	cout << "The array has now been sorted!\n\n";

	for (int i = 0; i < size; i++) { outfile << list[i] << endl; }

	cout << "The output has been written into the output file!\n";

	outfile << "Number of comparisons made using QuickSort: " << comparison_counter << endl;
	infile.close();
	outfile.close();

	infile.open("input.txt");
	if (!infile) { cout << "Error: The input file could not be opened" << endl; }
	else { cout << "The input file was successfully opened!" << endl; }

	comparison_counter = 0;

	infile >> size;
	cout << "There are " << size << " elements in the array" << endl;

	double *list2 = new double[size];

	for (int i = 0; i < size; i++) {
		infile >> list2[i];
	}

	ofstream outfile2;
	outfile2.open("rmk0017_merge.txt");

	
	cout << "**************************\n";
	cout << "Implementing MergeSort now...\n";
	MergeSort(list2, size);
	cout << "**************************\n";
	cout << "The array has now been sorted!\n";
	for (int i = 0; i < size; i++) { outfile2 << list2[i] << endl; }
	outfile2 << "Number of comparisons made using MergeSort : " << comparison_counter << endl;
	

	infile.close();
	outfile2.close();
}

bool LESS(double a, double b)
{
	comparison_counter++;
	if (a < b) { return true; }
	else { return false; }
}

void swap(double* a, double* b)
{
	double temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

void QuickSort(double A[], int left, int right, int n)
{
	/*
	if (left < right) {

		int pivot = Partition(A, left, right);
		QuickSort(A, left, pivot - 1, n); //Quicksorting the left side of the pivot
		QuickSort(A, pivot + 1, right, n); //Quicksorting the right side of the pivot

	}*/

	//My recursive implementation failed everytime the input was a large array and the worst case for quicksort as my computer would reach the limit of the recursion buffer
	//So I decided to employ an iterative implementation 
	
	stack<double> hold;
	hold.push(left);
	hold.push(right);

	while (!hold.empty()) {
		right = hold.top();
		hold.pop();
		left = hold.top();
		hold.pop();

		int pivot = Partition(A, left, right);

		if ((pivot - 1) > left){
			hold.push(left);
			hold.push(pivot - 1);
		}

		if ((pivot + 1) < right) {
			hold.push(pivot + 1);
			hold.push(right);
		}
	}
	
	
}

int Partition(double A[], int left, int right)
{
	int index = left;

	for (int i = left; i < right; i++) {
		if (LESS(A[i], A[right])) { swap(&A[i], &A[index]); index++; }
	}

	swap(&A[right], &A[index]);
	return index;
}

void MergeSort(double A[], int n)
{
	if (n > 1) {
		int m = ceil((double)n / 2);
		
		int q = n - m;

		double *B = new double[m];
		double *C = new double[q];

		for (int i = 0; i <= m; i++) {
			B[i] = A[i];
		}

		for (int i = 0, j = m; i <= (q); i++, j++) {
			C[i] = A[j];
		}

		
		MergeSort(B, m);
		MergeSort(C, q);
		Merge(B, C, A, m, q);	
	}
}

void Merge(double B[], double C[], double A[], int p, int q)
{
	int i = 0;
	int j = 0;

	for (int k = 0; k <= (p + q) - 1; k++) {
		if (i >= p) {
			A[k] = C[j]; j++;
		}
		else if (j >= q) {
			A[k] = B[i]; i++;
		}
		else {
			if (!(LESS(C[j], B[i]))) {
				A[k] = B[i]; i++;
			}
			else {
				A[k] = C[j]; j++;
			}
		}
	}
}

