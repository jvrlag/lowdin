// Computing fluctuations in particle number
// 160215-190325
#include"Lowdin.h"
#include"Manybody.h"
#include"Fermions.h"



int main(int argc, char *argv[])
{
     Rand_Open(0);
     long N=40; // size of initial matrix
     Input(N,"-N",argc,argv);
     long Ne=N/2;

     double alfa=0.5;
     Input(alfa,"-alpha",argc,argv);
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
	  printf("%5g %16.12g %16.12g\n",x,real(ntot),sqrt(real(n2-norm(ntot))));

	  

     }
     
}
