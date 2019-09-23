// 190307
// Lowdin classes

#include"Lowdin.h"

Box_Determ::Box_Determ()
{
     computed=false;
     Tol=1e-14;
}

void Box_Determ::Prepare(const MatrixC &a)
{
     Copy(A,a);
     N=A.N1;
     detA=Det(A);
     //printf("Det(A): %16.12g\n",abs(detA));
     if (abs(detA)>Tol)
     {
//	  printf("Matrix A is regular\n");
	  singular=false;
	  Ainv=Invert(A);
	  computed=true;
//	  printf("Done\n");
	  return;
     }
     singular=true;
     // The SVD trick upon Löwdin
     A.SVD(U,V,Sigma);
     //Sigma.Write();

     // Detect zeros
     Z=0;
     for (long i=N;i>=1;i--)
	  if (Sigma(i)<Tol) Z++;
	  else break;

     cmplx signA=detA/abs(detA);
     Theta=signA*Sigma.Prod(1,N-Z);
     
     Ainv.Create(N);
     Ainv.Zero();
     for (long i=1;i<=N-Z;i++)
	  Ainv+=(1./Sigma(i))*Ket_Bra(V.Row(i),U.Col(i));
     computed=true;
}


cmplx Box_Determ::Compute(const MatrixC &B, const MatrixC &C, const MatrixC &D) const
{
     if (!computed) Error("Not ready to compute!!");

     if (!singular)
     {
	  MatrixC Ared=D-C*Ainv*B;
	  return detA*Det(Ared);
     }
     
     MatrixC Dp=D-C*Ainv*B;
     MatrixC Ub=-Herm(B)*Part(U,1,N-Z+1,N,N);
     MatrixC Vb=C*Herm(Part(V,N-Z+1,1,N,N));

     //printf("Det(Dp): %16.12g\n",abs(Det(Dp)));
     if (abs(Det(Dp))<Tol) return 0.0;
     
     return Theta*Det(Dp)*Det(Herm(Ub)*Invert(Dp)*Vb);
}

Lowdin::Lowdin(const MatrixC &u, const MatrixC &v) : U(u), V(v)
{
     N=U.N1;
     Ne=U.N2;
     MatrixC A(Ne);
     for (long i=1;i<=Ne;i++)
	  for (long j=1;j<=Ne;j++)
	       A(i,j)=::Dot(U.Col(i),V.Col(j));
     BD.Prepare(A);
}

cmplx Lowdin::Dot() const
{
     return BD.detA;
}

cmplx Lowdin::Elem(long I, long J) const   // <U|c_I c^+_J|V>
{
     MatrixC B, C, D;
     B.Create(Ne,1);
     for (long i=1;i<=Ne;i++)
     	  B(i,1)=conj(U(I,i));
     C.Create(1,Ne);
     for (long i=1;i<=Ne;i++)
	  C(1,i)=V(I,i);
     D.Create(1);
     D(1,1)=1.0;
     return BD.Compute(B,C,D);
}

// <U| c_I1 c_I2 c+_J1 c+_J2 |V>
cmplx Lowdin::Elem(long I1, long I2, long J1, long J2) const
{
     MatrixC B, C, D;
     B.Create(Ne,2);
     for (long i=1;i<=Ne;i++)
     {
	  B(i,1)=conj(U(I1,i));
	  B(i,2)=conj(U(I2,i));
     }
     C.Create(2,N/2);
     for (long i=1;i<=N/2;i++)
     {
	  C(1,i)=V(J1,i);
	  C(2,i)=V(J2,i);
     }
     D.Create(2);
     D(1,1)=(I1==J1);
     D(1,2)=(I1==J2);
     D(2,1)=(I2==J1);
     D(2,2)=(I2==J2);
     return BD.Compute(B,C,D);     
}

cmplx Lowdin::NiNj(long I, long J) const
{
     // first term: <U|V>
     cmplx term0=BD.detA;
//     printf("Dot: %16.12g\n",real(term0));     

     // Now, using the routines
     MatrixC B, C, D;
     B.Create(N/2,1);
     for (long i=1;i<=N/2;i++)
     	  B(i,1)=conj(U(I,i));
     C.Create(1,N/2);
     for (long i=1;i<=N/2;i++)
	  C(1,i)=V(I,i);
     D.Create(1);
     D(1,1)=1.0;
     cmplx term1=BD.Compute(B,C,D);

     for (long i=1;i<=N/2;i++)
     	  B(i,1)=conj(U(J,i));
     for (long i=1;i<=N/2;i++)
	  C(1,i)=V(J,i);
     D(1,1)=1.0;
     cmplx term2=BD.Compute(B,C,D);

     B.Create(N/2,2);
     for (long i=1;i<=N/2;i++)
     {
	  B(i,1)=conj(U(I,i));
	  B(i,2)=conj(U(J,i));
     }
     C.Create(2,N/2);
     for (long i=1;i<=N/2;i++)
     {
	  C(1,i)=V(I,i);
	  C(2,i)=V(J,i);
     }
     D.Create(2);
     D(1,1)=1.0;
     D(1,2)=(I==J);
     D(2,1)=(I==J);
     D(2,2)=1.0;
     cmplx term3=BD.Compute(B,C,D);
     
     cmplx theor=term0-term1-term2+term3;
     if (I==J) theor+=term2;
     // printf("Theor: %16.12g %16.12g\n",real(theor),imag(theor));
     return theor;
}
     
