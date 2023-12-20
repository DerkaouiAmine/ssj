package WafomExperiments;


import java.util.Arrays;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;



/**
 * This class tests the frequency at which, if we perform a Left Matrix Scramble on our generating matrix,
 * we might end up with the same generator matrix. In such a case, an LMS (Left Matrix Scramble) is not necessary.
 * The tests here show that for \( w > k \geq 4 \) in any dimension, no matrices after an LMS are equal to the original\( C_j \).
 * This demonstrates the robustness of the method.
 */


public class NoNeedForLMS{
	public static void main(String[] args) 
	{
		int dim=3;
		int k=14;
		int w=31;
		int nbreSimulate=100000000;
		int[] [] IdentDim = new int [nbreSimulate][];


		for (int i=0;i<nbreSimulate;i++) {
			DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
			int[] a1=Sob.getGeneratorMatricesTrans();
			Sob.leftMatrixScramble(new MRG32k3a());
			int[] a2=Sob.getGeneratorMatricesTrans();
			IdentDim[i]=CompareMatrices(a1,a2,k,dim);
		}



		for (int i=0;i<IdentDim.length;i++)

		{

			for (int j=0;j<IdentDim[0].length;j++)
			{

				if (IdentDim[i][j]!=0)
				{
					System.out.println("No need for LMS in :");
					System.out.println("Simulation number :"+i+"\n");
					System.out.println("for Dimension number :"+j+"\n");

				}


			}}

		/*	for (int i=0;i<IdentDim.length;i++)

			{

for (int j=0;j<IdentDim[0].length;j++)
{System.out.print(IdentDim[i][j]+"  ");
	}

			System.out.println();}*/



		System.out.print(IdentDim.length+"  "+	 IdentDim[0].length);




	}

	private static int[] CompareMatrices(int[] tableau1 , int[] tableau2,int k,int dim) {

		//int k = 3; // Taille des sous-tableaux

		int l=0;
		//int numMatrice=0;
		int [] numDim = new int[dim];
		for (int i = 0; i < (tableau1.length)/dim; i +=k) {

			int[] subArray1 = Arrays.copyOfRange(tableau1, i*k, (i*k)+(k-1));
			int[] subArray2 = Arrays.copyOfRange(tableau2, i*k, (i*k)+(k-1) );

			if (Arrays.equals(subArray1, subArray2)) {
				numDim[l]=1;

			}
			else {
				numDim[l]=0; 
			}
			l+=1;
		}

		return numDim;
	}

}



















