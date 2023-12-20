package WafomExperiments;

import java.awt.datatransfer.SystemFlavorMap;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;





/**
 * In this class, we naively search for point sets with low WAFOM by generating several generating matrices and calculating the WAFOM each time. We retain only the point sets with the best WAFOMs.
 */

public class NaifLowWafomSearch {
	public static void calc(int k) { 
		long startTime = System.currentTimeMillis();
		System.out.println("je calcule pour k="+k);
		int dim=16;
		int dimension=0;
		//int k=14;//k=23 est le maximum dans notre cas
		int w=30;


		double functionValue=0;

		String chemin="/home/derkaoui/Desktop/ResultatsSob";
		RandomStream randomStream = new MRG32k3a();


		Tally stats=new Tally("Simul");
		Tally stats1=new Tally("Vari");

		int nbreSimulateLMS=1000;

		//int N=(int) Math.pow(2, k);//equivalent a 1<<k

		double [] WafomResult=new double[nbreSimulateLMS];

		int[][] []tableauPrincipalGenerLMS = new int[nbreSimulateLMS][][];

		double [][] pointSet=new double[1<<k][dim];



		for (int i=0;i<nbreSimulateLMS;i++) {
			stats.init();
			pointSet=new double[1<<k][dim];

			/**
			 * We generate our point sets 'nbreSimulateLMS' times.
			 */

			DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
			Sob.leftMatrixScramble(new MRG32k3a());
			/**
			 * We retain  generating matrices.
			 */

			tableauPrincipalGenerLMS[i] =  Sob.getGeneratorMatrices(k);
			/**
			 * We calculate the WAFOM each time.
			 */

			Wafom waf=new Wafom (Sob,1,dim,w,k);

			WafomResult[i]=waf.calcWafom();


		}


		/**
		 * We rearrange our results.
		 */


		int [] result1=findMinMax(WafomResult);
		double []rearrangeWafom = arrangeVecteur(WafomResult);

		System.out.println("Le plus petit Wafom : "+ WafomResult[result1[0]]);
		System.out.println("Le plus grand Wafom : "+ WafomResult[result1[1]]);

		System.out.println("Best gene Matri LMS");


		for (int i=0;i<tableauPrincipalGenerLMS[0].length;i++)
		{
			for (int j=0;j<tableauPrincipalGenerLMS[0][0].length;j++)
			{
				System.out.print(tableauPrincipalGenerLMS[result1[0]][i][j]+" ");
			}

			System.out.println();
		}

		System.out.println("Worst gene Matri LMS");


		for (int i=0;i<tableauPrincipalGenerLMS[0].length;i++)
		{
			for (int j=0;j<tableauPrincipalGenerLMS[0][0].length;j++)
			{
				System.out.print(tableauPrincipalGenerLMS[result1[1]][i][j]+" ");
			}

			System.out.println();

		}







		long endTime = System.currentTimeMillis();
		long executionTime = endTime - startTime;

		System.out.println("Temps d'exécution : " + executionTime + " millisecondes");

	}
	public static void main(String[] args) throws FileNotFoundException { 


		for (int k=1;k<23;k++) {
			calc(k);
		}
	}




	/**
	 * We find the minimum and the maximum.
	 */

	public static int[] findMinMax(double [] v)
	{
		int minIndex = 0;
		int maxIndex = 0;
		double minValue = v[0];
		double maxValue = v[0];

		for (int i = 1; i < v.length; i++) {
			if (v[i] < minValue) {
				minValue = v[i];
				minIndex = i;
			}

			if (v[i] > maxValue) {
				maxValue = v[i];
				maxIndex = i;
			}
		}
		int[] indices = {minIndex, maxIndex};
		return indices;
	}

	public static double[] arrangeVecteur(double[] vecteur) {
		double[] vecteurCopie = Arrays.copyOf(vecteur, vecteur.length);

		Arrays.sort(vecteurCopie);

		return vecteurCopie;
	}



}




