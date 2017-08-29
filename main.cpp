/*
 * main.c
 * 
 * Copyright 2017 Christian <wurzelblatt@ubuntu-vm>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <stdio.h>
#include "nr3.h"
#include "ludcmp.h"
#include <fstream>
#include <string.h>

using namespace std;

//declare and code function before it is called in main!
void print_mat(MatDoub M)
{
	Int n = M.nrows();
	Int m = M.ncols();
	Doub big, temp;
	MatInt Mdez(m,n);
	temp = 0.0;
	big = 0.0;
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			Mdez[i][j] = (int)(log10(abs(M[i][j]))+1);
			if (Mdez[i][j] <= -1) Mdez[i][j] += (int)(-log10(abs(M[i][j])));
			if (Mdez[i][j] <= 0 && Mdez[i][j] >= -1) Mdez[i][j] += 1;
			if(M[i][j] == 0) Mdez[i][j] = 1;
			if(M[i][j] - (int)M[i][j] != 0) Mdez[i][j] +=2;
			if(M[i][j] < 0) Mdez[i][j] +=1;
		}
	}
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			if((temp=abs(Mdez[i][j])) > big) big=temp;
		}
	}
	
	Int nDez = big;
	for (int i=0; i<n; i++) {
		cout<<"|  ";
		for (int j=0; j<m; j++) {
			char CharArr[nDez];
			int it = (nDez)-Mdez[i][j];
			for (int k=0; k<it; k++) {
				CharArr[k] = ' ';
			}
			
			int itemp = 0;
			if (M[i][j] - (int)M[i][j] == 0) {		
				for (int k=nDez-1; k>=it; k--) {
					int mod = pow(10, (nDez-k));
					itemp = ((int)abs(M[i][j]) - itemp) % mod;
					CharArr[k] = ((itemp/(mod/10))+48);
				}
			}
			else {
				for (int k=nDez-3; k>=it; k--) {
					int mod = pow(10, (nDez-k-2));
					itemp = ((int)abs(M[i][j]) - itemp) % mod;
					CharArr[k] = ((itemp/(mod/10))+48);
					CharArr[nDez-2] = '.';
					int idig = int((M[i][j] - int(M[i][j])) * 10);
					CharArr[nDez-1] = abs(idig)+48;
				}
			}
			
			if(M[i][j] < 0) CharArr[it] = '-';
			//else CharArr[it-1] = ' ';
			
			for (int k = 0; k<nDez; k++) {
				cout<<CharArr[k];
			}
			cout<<"  ";	
		}
		cout<<'|';
		cout<<endl;		
	}
}

void print_vec(VecDoub V)
{	
	Int n = V.size();
	Doub temp, big;
	temp = 0.0;
	big = 0.0;
	
	//VecInt Vdez(n);
	int Vdez[n];
	for (int i=0; i<n; i++) {
		Vdez[i] = (int)(log10(abs(V[i]))+1);
		if (Vdez[i] <= -1) Vdez[i] += (int)(-log10(abs(V[i])));
		if (Vdez[i] <= 0 && Vdez[i] >= -1) Vdez[i] += 1;
		if(V[i] == 0) Vdez[i] = 1;
		if(V[i] - (int)V[i] != 0) Vdez[i] += 2;
		if (V[i] < 0) Vdez[i] += 1;
	}
	
	for (int i=0; i<n; i++) {
		if((temp=abs(Vdez[i])) > big) big=temp;
	}
	
	Int nDez = big;
	for (int i=0; i<n; i++) {
		cout<<"|  ";
		char CharArr[nDez];
		int it = (nDez)-Vdez[i];
		for (int k=0; k<it; k++) {
			CharArr[k] = ' ';
		}
		
		
		int itemp = 0;	
		if (V[i] - int(V[i]) == 0) {
			for (int k=nDez-1; k>=it; k--) {
				int mod = pow(10, (nDez-k));
				itemp = ((int)abs(V[i]) - itemp) % mod;
				CharArr[k] = ((itemp/(mod/10))+48);
			}
		}
		else {
			for (int k=nDez-3; k>=it; k--) {
				int mod = pow(10, (nDez-k-2));
				itemp = ((int)abs(V[i]) - itemp) % mod;
				CharArr[k] = ((itemp/(mod/10))+48);
				CharArr[nDez-2] = '.';
				int idig = int((V[i] - int(V[i])) * 10);
				CharArr[nDez-1] = abs(idig)+48;
			}						 
		}
		
		if(V[i] < 0) CharArr[it] = '-';
		//else CharArr[it-3] = ' ';
		
		for (int k = 0; k<nDez; k++) {
			cout<<CharArr[k];		
		}
		cout<<"  |";
		cout<<endl;	
	}


}

void print_vec(VecInt V)
{
	Int n= V.size();
	Doub temp, big;
	temp = 0.0;
	big = 0.0;
	int Vdez[n];
	for (int i=0; i<n; i++) {
		Vdez[i] = (int)(log10(abs(V[i]))+1);
		if(V[i] == 0) Vdez[i] = 1;
		if (V[i] < 0) Vdez[i] += 1;
	}
	for (int i=0; i<n; i++) {
		if((temp=abs(Vdez[i])) > big) big=temp;
	}
	Int nDez = big;
	for (int i=0; i<n; i++) {
		cout<<"|  ";
		char CharArr[nDez];
		int it = (nDez)-Vdez[i];
		for (int k=0; k<it; k++) {
			CharArr[k] = ' ';
		}
		int itemp = 0;	
		for (int k=nDez-1; k>=it; k--) {
			int mod = pow(10, (nDez-k));
			itemp = ((int)abs(V[i]) - itemp) % mod;
			CharArr[k] = ((itemp/(mod/10))+48);
		}
	if(V[i] < 0) CharArr[it] = '-';
	for (int k = 0; k<nDez; k++) {
		cout<<CharArr[k];		
	}
	cout<<"  |";
	cout<<endl;
	}
}	
		



int main()
{
	//MatDoub A(4,4);
	//Doub elements = 0.;
	//for (int i = 0; i < 4; i++) {
		//for (int j = 0; j < 4; j++) {
			//elements = elements + 1.;
			//A[i][j] = elements;
		//}
	//}
	const Int n = 3;
	
	MatDoub AA(n,n);
	AA[0][0] = 1;
	AA[0][1] = 2;
	AA[0][2] = 4;
	AA[1][0] = 3;
	AA[1][1] = 8;
	AA[1][2] = 14;
	AA[2][0] = 2;
	AA[2][1] = 6;
	AA[2][2] = 13;

	
	//MatDoub AA(n,n);
	//AA[0][0] = 2;
	//AA[0][1] = -6;
	//AA[0][2] = -1;
	//AA[0][3] = -1;
	//AA[1][0] = -3;
	//AA[1][1] = -1;
	//AA[1][2] = 7;
	//AA[1][3] = 5;
	//AA[2][0] = -8;
	//AA[2][1] = 1;
	//AA[2][2] = -2;
	//AA[2][3] = 15;
	//AA[3][0] = 1;
	//AA[3][1] = 10;
	//AA[3][2] = 10;
	//AA[3][3] = -4;
	
	//VecDoub V(n);
	//V[0] = 0.023;
	//V[1] = 0.375;
	//V[2] = -14;
	//V[3] = 0.2391;
	
	VecDoub b(n);
	b[0] = 3;
	b[1] = 13;
	b[2] = 4;
	//b[3] = 5;
	
	VecDoub x(n);
	
	LUdcmp my_ludcmp(AA);
	
	my_ludcmp.solve(b,x);
		
	
	//cin.get();
	cout<<"Matrix A:"<<endl;
	print_mat(AA);		 
	cin.get();
	

	//cout<<"Vector b:"<<endl;
	//print_vec(b);
	//cin.get(); 

	//cout<<"Size of Indx: "<<my_ludcmp.indx.size()<<endl; //depends on how often we interchanged rows
	cout<<"Indx: "<<endl;
	print_vec(my_ludcmp.indx);
	cin.get();
	
	
	//cout<<(my_ludcmp.lu[0][0])<<endl;
	//cout<<(my_ludcmp.lu[0][1])<<endl;
	//cout<<(my_ludcmp.lu[0][2])<<endl;
	//cout<<(my_ludcmp.lu[1][0])<<endl;
	//cout<<(my_ludcmp.lu[1][1])<<endl;
	//cout<<(my_ludcmp.lu[1][2])<<endl;
	//cout<<(my_ludcmp.lu[2][0])<<endl;
	//cout<<(my_ludcmp.lu[2][1])<<endl;
	//cout<<(my_ludcmp.lu[2][2])<<endl;
	//cin.get();
	cout<<(my_ludcmp.vv[0])<<endl;
	cout<<(my_ludcmp.vv[1])<<endl;
	cout<<(my_ludcmp.vv[2])<<endl;
	cout<<(x[0])<<endl;
	cout<<(x[1])<<endl;
	cout<<(x[2])<<endl;		
	cin.get();	

	
	cout<<"LU decomposition:"<<endl;
	print_mat(my_ludcmp.lu);
	cout<<endl;
	cout<<"Vector x:"<<endl;
	print_vec(x);
	cout<<endl;
	cout<<"Vector P:"<<endl;
	print_vec(my_ludcmp.P);

	cin.get();
	
	//system("pause");   //nur in Windows
	return 0;
}



