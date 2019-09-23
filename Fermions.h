// Functions for dealing with free-fermionic systems
// 150413-170519
#ifndef FREE_FERMIONS_H
#define FREE_FERMIONS_H
#include"Manybody.h"

double H2(double x);

// Hopping matrix with hoppings T, 1D system with PBC
MatrixC Hopping_Matrix(const Vector &T);

MatrixC Hopping_Matrix(const Graph &G, const Vector &T);

MatrixC FF_Corr(const MatrixC &U, long Ne);

// Returns the entanglement orbitals and energies for a certain block
void FF_Block_Entanglement(MatrixC &B, Vector &E, 
			   const MatrixC &Corr, const List &P);

double FF_Block_Entropy(const MatrixC &Corr, const List &P);

MatrixC FF_Dynamics(const MatrixC U0, 
		    const MatrixC U, const Vector E, double t);

cmplx FF_Dot(const MatrixC &U1, const MatrixC &U2, long Ne);

#endif


