// Computing fluctuations in particle number
// 160215-190325
#include"Lowdin.h"
#include"Manybody.h"
#include"Fermions.h"

// Compute the full wf of the Slater determinant
VectorC Slater_Det_2(const MatrixC &B, const List &K)
{
     long L=B.N1;
     long Nt=Pow_2(L);
     VectorC Psi(Nt);
     Psi(1)=1.0;
     for (long k=1;k<=K.N;k++)
     {
	  MatrixC Creator(Nt);
	  for (long i=1;i<=L;i++)
	       Creator+=B(i,K(k))*Ferm_Site_Op(C_Op(+1),i,L);
	  Psi=Creator*Psi;
     }
     return Psi;
}

double Occupation(const VectorC &Psi, const List &Lr, long N)
{
     MatrixC Op(Psi.N);
     for (long i=1;i<=Lr.N;i++)
	  Op+=Site_Op(C_Op(0),Lr(i),N);
     return real(Op.Elem(Psi,Psi));
}

double Occ_Sq(const VectorC &Psi, const List &Lr, long N)
{
     MatrixC Op(Psi.N);
     for (long i=1;i<=Lr.N;i++)
	  Op+=Site_Op(C_Op(0),Lr(i),N);
     Op=Op*Op;
     return real(Op.Elem(Psi,Psi));
}




int main()
{
     Rand_Open(0);
     long N=8; // size of initial matrix
     long Ne=N/2;

     double alfa=0.9;
     Vector J(N);
     J(N/2)=-1.0;
     for (long i=1;i<N/2;i++)
	  J(N/2-i)=J(N/2+i)=-pow(alfa,2*i-1);

     MatrixC Hop=Hopping_Matrix(J);
     MatrixC Basis; Vector Eigen;
     Hop.Diagonalize(Basis,Eigen);
     printf("# Energy: %16.12g\n",Eigen.Sum(1,Ne));
     printf("# Energy2: %16.12g\n",Eigen.Sum(1,Ne-1)+Eigen(Ne+1));
     
     MatrixC U=Part(Basis,1,1,N,Ne);
     MatrixC V=Part(Basis,1,1,N,Ne+1);
     V.Remove_Col(Ne);



     
     VectorC PsiU=Slater_Det_2(U,List_Range(1,Ne));
     VectorC PsiV=Slater_Det_2(V,List_Range(1,Ne));
     cmplx dot=Dot(PsiU,PsiV);
     printf("# Dot: %16.12g\n",real(dot));


     Lowdin L_UU(U,U), L_VV(V,V), L_UV(U,V);

     printf("# X Lowdin_occ Lowdin_sigma Ex_occ Ex_sigma Ex_entropy\n");
     
     for (double x=0.0;x<=1.0;x+=0.01)
     {
	  double alpha=sqrt(1.-x), beta=sqrt(x);
     
	  cmplx ntot=0.0;
	  for (long i=1;i<=N/2;i++)
	  {
	       ntot+=norm(alpha)*(1.-L_UU.Elem(i,i));
	       ntot+=norm(beta)*(1.-L_VV.Elem(i,i));
	       ntot+=2*real(conj(alpha)*beta*(L_UV.Dot()-L_UV.Elem(i,i)));
	  }

	  cmplx n2=0.0;
	  for (long i=1;i<=N/2;i++)
	       for (long j=1;j<=N/2;j++)
	       {
		    n2+=norm(alpha)*L_UU.NiNj(i,j);
		    n2+=norm(beta)*L_VV.NiNj(i,j);
		    n2+=2*real(conj(alpha)*beta*L_UV.NiNj(i,j));
	       }
	  printf("%5g %16.12g %16.12g ",x,real(ntot),sqrt(real(n2-norm(ntot))));

	  VectorC Psi=alpha*PsiU+beta*PsiV;
	  double occ=Occupation(Psi,List_Range(1,N/2),N);
	  double occ2=Occ_Sq(Psi,List_Range(1,N/2),N);

	  MatrixC Rho=Ket_Bra(Psi,Psi);
	  MatrixC RhoA=Trace_On(Rho,List_Range(1,N/2));
	  double S=Von_Neumann(RhoA);
	  
	  printf("%16.12g %16.12g %16.12g\n",occ,sqrt(occ2-Sqr(occ)),S);
//	  printf("Occ2: %16.12g\n",occ2);
	  

     }
     
}
