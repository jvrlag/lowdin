// Lowdin classes
#ifndef LOWDIN_H
#define LOWDIN_H

#include"MatrixC.h"

class Box_Determ
{
public:
     MatrixC A;
     bool computed;
     bool singular;
     MatrixC U, V; // SVD decomposition
     Vector Sigma;
     cmplx Theta; 
     long N; // size of initial system
     long Z; // number of zeros in singular values

     double Tol=1e-14; // tolerance to detect zeros
     cmplx detA; // determinant of A
     MatrixC Ainv; // inverse or pseudoinverse, depending!

     Box_Determ();
     ~Box_Determ() {};
     void Prepare(const MatrixC &a);
     cmplx Compute(const MatrixC &b, const MatrixC &c, const MatrixC &d) const;
};

class Lowdin
{
public:
     MatrixC U, V; // |U> and |V>
     long N; // size of the system
     long Ne; // number of particles
     Box_Determ BD;

     Lowdin(const MatrixC &u, const MatrixC &v);
     ~Lowdin() {};
     cmplx Dot() const;
     cmplx Elem(long I, long J) const;
     cmplx Elem(long I1, long I2, long J1, long J2) const;
     cmplx NiNj(long I, long J) const;
};
     
	  

#endif
