package wafomExperiments;

import java.awt.datatransfer.SystemFlavorMap;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

import umontreal.ssj.discrepancy.Wafom;

import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.discrepancy.WafomRQMC;
import umontreal.ssj.discrepancy.WafomRQMCFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.NiedXingSequenceBase2;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

/**
 * This class focuses on identifying point sets with low WAFOM. It achieves this
 * by applying a Left matrix scramble several times to each generating matrix,
 * and subsequently calculating the WAFOM. Only the point sets with the best
 * WAFOM values are retained. this Algorithme was described in My pdf and
 * in @cite{Shin Harase. Quasi-Monte Carlo point sets with small t-values and
 * WAFOM. Applied Mathematics and Computation, 254:318–326, 2015}
 */

public class LowWafomSearchAllGeneratingMatrixAtOnce {
	public static void calc(int k) {
		long startTime = System.currentTimeMillis();
		System.out.println("je calcule pour k=" + k);
		int dim = 12;
		int dimension = 0;
		int w = 30;
		int nbreSimulateLMS = 10000;
		int[][] tablGenerMatrices = new int[nbreSimulateLMS][];
		double functionValue = 0;

		String chemin = "/home/derkaoui/Desktop/ResultatsSob";
		RandomStream randomStream = new MRG32k3a();

		Tally stats = new Tally("Simul");
		Tally stats1 = new Tally("Vari");
		double[] WafomResult = new double[nbreSimulateLMS];

		int[][] tableauPrincipalGenerLMS = new int[nbreSimulateLMS][];
		int[][][] tableauPrincipalLeftMatrices = new int[nbreSimulateLMS][][];
		for (int i = 0; i < nbreSimulateLMS; i++) {
			stats.init();

			/**
			 * We generate our point sets 'nbreSimulateLMS' times.
			 */

			DigitalNetBase2 Sob = new SobolSequence(k, w, dim);
			// DigitalNetBase2 NX= new NiedXingSequenceBase2(k,w,dim);

			// we store the column "d" from L
			tableauPrincipalLeftMatrices[i] = Sob.leftMatrixScrambleBis(new MRG32k3a());

			WafomFast waf = new WafomFast(Sob, w, 3);

			// For Goda FOM
			// WafomRQMCFast waf1=new WafomRQMCFast (Sob,k,w,dim,3,new RandomShift(new
			// MRG32k3a()));

			tableauPrincipalGenerLMS[i] = Sob.getGeneratorMatricesTrans();

			WafomResult[i] = waf.calcWafom();

		}

		/**
		 * We rearrange our results.
		 */

		int[] result1 = findMinMax(WafomResult);
		double[] rearrangeWafom = arrangeVecteur(WafomResult);

		System.out.println("Le plus petit Wafom : " + WafomResult[result1[0]]);
		System.out.println("Le plus grand Wafom : " + WafomResult[result1[1]]);

		System.out.println("Best L in LMS");

		for (int i = 0; i < tableauPrincipalLeftMatrices[0].length; i++) {
			for (int j = 0; j < tableauPrincipalLeftMatrices[0][0].length; j++) {
				System.out.print(tableauPrincipalLeftMatrices[result1[0]][i][j] + " ");
			}

			System.out.println();
		}

		System.out.println("Worst L in LMS");

		for (int i = 0; i < tableauPrincipalLeftMatrices[0].length; i++) {
			for (int j = 0; j < tableauPrincipalLeftMatrices[0][0].length; j++) {
				System.out.print(tableauPrincipalLeftMatrices[result1[1]][i][j] + " ");
			}

			System.out.println();

		}

		System.out.println("generatring matrix \tilde{C}:");

		System.out.println("Best ");

		System.out.println(Arrays.toString(tableauPrincipalGenerLMS[result1[0]]));

		System.out.println("Worst");

		System.out.println(Arrays.toString(tableauPrincipalGenerLMS[result1[1]]));

		long endTime = System.currentTimeMillis();
		long executionTime = endTime - startTime;

		System.out.println("Temps d'exécution : " + executionTime + " millisecondes");

	}

	public static void main(String[] args) throws FileNotFoundException {

		for (int k = 1; k < 21; k++) {
			calc(k);
		}
	}

	/**
	 * We find the minimum and the maximum.
	 */

	public static int[] findMinMax(double[] v) {
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
		int[] indices = { minIndex, maxIndex };
		return indices;
	}

	public static double[] arrangeVecteur(double[] vecteur) {
		double[] vecteurCopie = Arrays.copyOf(vecteur, vecteur.length);

		Arrays.sort(vecteurCopie);

		return vecteurCopie;
	}

}
