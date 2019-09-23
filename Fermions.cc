// 150413
// Routines for free-fermion systems
#include"Fermions.h"

double H2(double x)
{
     if ((x<1e-30) || ((1.0-x)<1e-30))
	  return 0.0;
     return -(x*log(x)+(1.0-x)*log(1.0-x));
}

// Hopping matrix with hoppings T, 1D system with PBC
MatrixC Hopping_Matrix(const Vector &T)
{
     long L=T.N;
     MatrixC Hop(L);
     for (long k=1;k<L;k++)
          Hop(k,k+1)=Hop(k+1,k)=-T(k);
     Hop(1,L)=Hop(L,1)=-T(L);
     return Hop;
}

// Hopping matrix on a graph
MatrixC Hopping_Matrix(const Graph &G, const Vector &J)
{
     long L=G.N;
     MatrixC Hop(L);
     for (long k=1;k<=G.Nl;k++)
     {
	  long s1, s2;
	  G.Link_Sites(s1,s2,k);
	  Hop(s1,s2)=Hop(s2,s1)=-J(k);
     }
     return Hop;
}

MatrixC FF_Corr(const MatrixC &U, long Ne)
{
     long N=U.N1;
     MatrixC C(N);
     for (long k=1;k<=Ne;k++)
	  for (long i=1;i<=N;i++)
	       for (long j=1;j<=N;j++)
		    C(i,j)+=conj(U(i,k))*U(j,k);
     return C;
}

void FF_Block_Entanglement(MatrixC &B, Vector &E, 
			   const MatrixC &Corr, const List &P)

{
     long n=P.N;
     MatrixC Csub(n);
     for (long i=1;i<=n;i++)
     {
	  long i1=P(i);
	  Csub(i,i)=Corr(i1,i1);
	  for (long j=i+1;j<=n;j++)
	       Csub(i,j)=Csub(j,i)=Corr(i1,P(j));
     }
     Csub.Diagonalize(B,E);
}

double FF_Block_Entropy(const MatrixC &Corr, const List &P)
{
     MatrixC B; Vector Eigen;
     FF_Block_Entanglement(B,Eigen,Corr,P);
     double S=0.0;
     for (long i=1;i<=Eigen.N;i++)
	  S+=H2(Eigen(i));
     return S;
}

MatrixC FF_Dynamics(const MatrixC U0, 
		    const MatrixC U, const Vector E, double t)
{
     long L=U.N1;
     MatrixC Exp(L);
     for (long i=1;i<=L;i++)
	  Exp(i,i)=exp(-M_I * E(i) * t);
     return U * Exp * Herm(U) * U0;
}

// Compute <U1,Ne|U2,Ne>
cmplx FF_Dot(const MatrixC &U1, const MatrixC &U2, long Ne)
{
     MatrixC U=Herm(U1)*U2;     
     MatrixC Beta(Ne,Ne);
     for (long i=1;i<=Ne;i++)
      	  for (long j=1;j<=Ne;j++)
      	       Beta(i,j)=U(i,j);
     return Det(Beta);
}
