package wafomExperiments;

import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;

/**
 * This algorithm for searching low WAFOM was presented by Shin Harase. In:
 * Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357. It presents an
 * efficient method to obtain a low-WAFOM. It proceeds gradually and extends the
 * point set \(P_d = \{\mathbf{x}_0, \dots, \mathbf{x}_{2^d-1}\} \supsetneq
 * P_{d-1}\) for \(0 \leq d \leq k\), where each \(P_d\) is a subset of \(P\)
 * and \(|P|=2^k\). The algorithm tests for each \(d\) a large number of
 * different point sets, keeping the best ones in terms of WAFOM.
 * 
 * We operate on the matrix \tilde{C} = L * C, where L represents the left
 * matrix scramble. Randomizations are incrementally applied to each column of L
 * as we construct them.
 * 
 */

public class ExtendPointSetLowWafomBestColumn {// je prend par colonne
	public static void main(String[] args) throws FileNotFoundException {
		long startTime = System.currentTimeMillis();
		int dim = 5;
		int k = 21;
		int w = 30; // We set w=30 to be able to use the WAFOM Fast methods with w=30 and q=3.

		/**
		 * This refers to the number of randomizations assigned to each column of the
		 * generator matrix.
		 */
		int nbreTest = 10000;
		int[][] tableauPrincipalGenerLMS = new int[2 * k][];
		double[] bestWAFOMS = new double[k];
		double[] worstWAFOMS = new double[k];
		int[][] bestColonnes = new int[k][dim];
		int[][] worstColonnes = new int[k][dim];
		/**
		 * This is the generator matrix that we will modify as we progress with our
		 * tests.
		 */
		int[][] BEstGeneramatrix = new int[k][dim];
		int[] SobOriginale = new int[dim * k];
		int[][] OptimalLeftMatrix = new int[dim][w];
		// We take point sets from 2 points to 2^23 points.
		for (int d = 0; d < k; d++) {
			System.out.println("Je calcule pour k=" + (d + 1));
			double smallestWafom = (1 << k);
			double biggestWafom = -1;
			int[] colonne1 = new int[dim];
			int[] ColonneGenrat = new int[dim];

			for (int i1 = 0; i1 < nbreTest; i1++) {
				DigitalNetBase2 Sob = new SobolSequence(d + 1, w, dim);

				SobOriginale = Sob.getGeneratorMatricesTrans();
				// Here, we test the "d"-th column of the left triangular matrix "nbreTest"
				// times.

				/*
				 * if (d==0 && i1==0) { OptimalLeftMatrix=Sob.leftMatrixScrambleTrue1(new
				 * MRG32k3a())}; // This won't affect the result because, for the initial
				 * iteration, the left matrix scramble is fixed to all ones.
				 */

				OptimalLeftMatrix = Sob.leftMatrixScrambleExtendBis(new MRG32k3a(), OptimalLeftMatrix, d);

				// Extracting the d-th column of matrix L.

				colonne1 = extraireColonne(OptimalLeftMatrix, d);// ici i1 a la place de 0
				ColonneGenrat = Sob.getGeneratorMatricesTrans();
				/**
				 * We calculate the WAFOM each time
				 */
				// RandomShift rand = new RandomShift(new MRG32k3a());

				// WafomRQMCAccel lowWaf=new WafomRQMCAccel (Sob,1,dim,w,d+1,3,rand);

				WafomFast lowWaf = new WafomFast(Sob, w, 3);
				double currentWafom = lowWaf.calcWafom();

				/**
				 * We track the best and worst WAFOM values along with their corresponding best
				 * and worst columns of matrix L, and the generating matrices.
				 */

				if (currentWafom < smallestWafom) {
					smallestWafom = currentWafom;
					bestWAFOMS[d] = smallestWafom;
					bestColonnes[d] = colonne1;
					BEstGeneramatrix[d] = ColonneGenrat;
				}
				if (currentWafom > biggestWafom) {
					biggestWafom = currentWafom;
					worstWAFOMS[d] = biggestWafom;
					worstColonnes[d] = colonne1;
				}
			}
			long endTime = System.currentTimeMillis();
			long executionTime = endTime - startTime;

			System.out.println(executionTime);

			/**
			 * We keep the left matrix with the best column that generates the best points
			 * in terms of WAFOM to extend our set to a larger d.
			 */
			OptimalLeftMatrix = replaceColumn(OptimalLeftMatrix, bestColonnes[d], d);

			System.out.println("best WAFOMS");
			System.out.println(Arrays.toString(bestWAFOMS));

			System.out.println("BEST COLOMNS L");
			for (int i = 0; i < bestColonnes[0].length; i++) {
				for (int j = 0; j < bestColonnes.length; j++) {
					System.out.print(bestColonnes[j][i] + " ");
				}
				System.out.println();
			}

			System.out.println("worst WAFOMS L");
			System.out.println(Arrays.toString(worstWAFOMS));

			System.out.println("worst COLOMNS");
			for (int i = 0; i < worstColonnes[0].length; i++) {
				for (int j = 0; j < worstColonnes.length; j++) {
					System.out.print(worstColonnes[j][i] + " ");
				}
				System.out.println();
			}
			System.out.println("Best tilde Generatrice L*C:");
			System.out.println(Arrays.toString(BEstGeneramatrix[d]));

		}
	}

	/**
	 * This method is essential because we generate generator matrices with k
	 * columns, and we extract just the k-th column 'nbreTest' times to test it.
	 */

	public static int[] extraireColonne(int[][] tableau, int j) {
		int nbLignes = tableau.length;
		int[] colonne = new int[nbLignes];

		for (int i = 0; i < nbLignes; i++) {
			colonne[i] = tableau[i][j];
		}

		return colonne;
	}

	/**
	 * When we extract our 'nbre sets' of columns, we should replace them in our
	 * generator matrix to test them in terms of WAFOM.
	 */

	public static int[][] replaceColumn(int[][] matrix, int[] newColumn, int columnIndex) {
		if (columnIndex < 0 || columnIndex >= matrix[0].length) {
			throw new IllegalArgumentException("L'index de colonne spécifié est invalide.");
		}

		if (newColumn.length != matrix.length) {
			throw new IllegalArgumentException(
					"Le nouveau vecteur colonne doit avoir la même taille que le nombre de lignes de la matrice.");
		}

		for (int i = 0; i < matrix.length; i++) {
			matrix[i][columnIndex] = newColumn[i];
		}

		return matrix;
	}

}