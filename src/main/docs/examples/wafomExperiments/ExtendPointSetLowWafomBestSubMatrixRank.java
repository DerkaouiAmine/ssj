package wafomExperiments;


import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.discrepancy.WafomFast;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;

/**
 * This algorithm for searching low WAFOM was presented by Shin Harase.
 * In: Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357.
 * It presents an efficient method to obtain a low-WAFOM. It proceeds gradually
 * and extends the point set \(P_d = \{\mathbf{x}_0, \dots, \mathbf{x}_{2^d-1}\} \supsetneq P_{d-1}\) for \(0 \leq d \leq k\),
 * where each \(P_d\) is a subset of \(P\) and \(|P|=2^k\). The algorithm tests for each \(d\) a large number of different point sets,
 * keeping the best ones in terms of WAFOM.
*/
/**
 * This class is similar to ExtendPointSetLowWafomBestColumn, but here our focus
 * is directly on constructing the columns of the generating matrix C bit by bit
 * at random.
 * 
 * We will work on the generator matrices, so all randomizations are done on the
 * columns of the generator matrices gradually for each d. We ensure that each
 * "d-1" x "d-1" submatrix is invertible using Gaussian Elimination. Otherwise,
 * since we work in \( \mathbf{F}_2 \), it suffices to swap the (d,d) element of
 * the new matrix.
 */

public class ExtendPointSetLowWafomBestSubMatrixRank {// je prend par colonne
	public static void main(String[] args) throws FileNotFoundException {
		long startTime = System.currentTimeMillis();

		int dim = 5;
		int k = 20;
		int w = 30;

		/**
		 * This refers to the number of randomizations assigned to each column of the
		 * generator matrix.
		 */
		int nbreTest = 10000;

		double[] bestWAFOMS = new double[k];
		double[] worstWAFOMS = new double[k];
		int[][] bestColonnes = new int[k][w * dim];
		int[][] worstColonnes = new int[k][dim];
		int[] bestcolumn = new int[w * dim];
		int[] BestGeneratorMatrices = new int[dim * k];
		int[][] BEstGeneramatrix = new int[k][dim];
		int[][] Generamatrix = new int[dim][0];
		int[][][] optimaleMatrix = new int[dim][w][k];
		for (int d = 0; d < k; d++) {
			System.out.println("Je calcule pour k=" + (d + 1));
			double smallestWafom = (1 << k);
			double biggestWafom = -1;

			// For each k, we test 'nbreTest' times.
			for (int i1 = 0; i1 < nbreTest; i1++) {
				if (d == 0) {
					optimaleMatrix = generateMatrix(dim, w, d + 1, 2, new MRG32k3a());

				} else {
					optimaleMatrix = generateMatrix1(dim, w, d + 1, 2, optimaleMatrix, new MRG32k3a());
				}

				/**
				 * Check the invertibility of submatrices within the optimaleMatrix.
				 * 
				 * - `allInvertible`: Flag indicating if all submatrices are invertible. -
				 * `problematicIndices`: Array storing indices of non-invertible submatrices. -
				 * `count`: Counter for the number of non-invertible submatrices.
				 */

				boolean allInvertible = true;
				int[] problematicIndices = new int[dim]; // stock the index of the non invertible sub-matrix .
				int count = 0;

				for (int j = 0; j < dim; j++) {
					if (!isSubMatrixInvertible1(optimaleMatrix[j], d + 1)) {
						allInvertible = false;
						problematicIndices[count++] = j;
					}
				}

				if (!allInvertible) {
					// fix the non invertible sub matrix by a swap of the element (d,d).
					for (int i = 0; i < count; i++) {
						int index = problematicIndices[i];
						if (optimaleMatrix[index][d][d] == 0) {
							optimaleMatrix[index][d][d] = 1;
						} else {
							optimaleMatrix[index][d][d] = 0;
						}
					}
				}

				// stock the new columns
				int index = 0;
				for (int j = 0; j < dim; j++) {
					int[] colonne1 = extraireColonne(optimaleMatrix[j], d);
					for (int i = 0; i < colonne1.length; i++) {
						bestcolumn[i + index] = colonne1[i];

					}
					index += w;

				}
				// and the generating matrices
				int[] generate = generatorMatricesFromStandardFormat(optimaleMatrix);

				DigitalNetBase2 dn = new DigitalNetBase2(d + 1, w, dim, generate);

				/**
				 * We calculate the WAFOM each time
				 */

				WafomFast lowWaf = new WafomFast(dn, w, 3);

				double currentWafom = lowWaf.calcWafom();

				/**
				 * We keep track of the best and worst WAFOM values.
				 */

				if (currentWafom < smallestWafom) {
					smallestWafom = currentWafom;

					bestWAFOMS[d] = smallestWafom;
					bestColonnes[d] = bestcolumn;
					BEstGeneramatrix[d] = generate;

				}

				if (currentWafom > biggestWafom) {
					biggestWafom = currentWafom;

					worstWAFOMS[d] = biggestWafom;

				}

			}
			long endTime = System.currentTimeMillis();
			long executionTime = endTime - startTime;

			System.out.println(executionTime);

			/**
			 * We keep the generator matrix with the best column that generates the best
			 * points in terms of WAFOM to extend our set to a larger d.
			 */
			int indexx = 0;
			for (int j = 0; j < dim; j++) {
				int[] Column = new int[w];
				for (int i = 0; i < w; i++) {

					Column[i] = bestColonnes[d][i + indexx];
				}
				indexx += w;

				optimaleMatrix[j] = replaceColumn(optimaleMatrix[j], Column, d);
			}
			System.out.println("best WAFOMS");
			System.out.println(Arrays.toString(bestWAFOMS));

			System.out.println("worst WAFOMS");
			System.out.println(Arrays.toString(worstWAFOMS));

			System.out.println("Best GeneratriceNormalement sa reste la meme et sa evolue");
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

	// Here, we generate the generating matrix bit by bit using randomization.

	public static int[][][] generateMatrix(int s, int w, int k, int b, RandomStream stream) {
		int[][][] L = new int[s][w][k];

		for (int j = 0; j < s; j++) {
			for (int l = 0; l < w; l++) {
				for (int c = 0; c < k; c++) {
					// LOI de PARETO

					// L[j][l][c] = stream.nextInt(0, b - 1);
					L[j][l][c] = (stream.nextDouble() < 0.8) ? 1 : 0; // 80% de chance de 1 et 20% de 0

				}
			}
		}
		return L;
	}

	// Here, we fix the previous "d-1" columns and generate the d-th column
	// randomly.

	public static int[][][] generateMatrix1(int s, int w, int k, int b, int[][][] r, RandomStream stream) {// r a k-1
																											// colonne
																											// de
																											// 0,...,k-2
		int[][][] L = new int[s][w][k];

		for (int j = 0; j < r.length; j++) {
			for (int l = 0; l < r[0].length; l++) {
				for (int c = 0; c < r[0][0].length; c++) {
					L[j][l][c] = r[j][l][c];
				}
			}
		}

		for (int j = 0; j < s; j++) {
			for (int l = 0; l < w; l++) {
				// L[j][l][k-1] = stream.nextInt(0, b - 1);

				L[j][l][k - 1] = (stream.nextDouble() < 0.8) ? 1 : 0;
			}
		}

		return L;
	}

	//////////////////////////////////////////////////////////////////////////

//Gauss elimination Method

	public static boolean isSubMatrixInvertible1(int[][] matrix, int k) {

		k = matrix[0].length;
		int[][] subMatrix = new int[k][k];

		// get the submatrix k*k
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {
				subMatrix[i][j] = matrix[i][j];
			}
		}

		// Apply Gauss
		for (int i = 0; i < k; i++) {
			// find the pivot
			int pivotRow = -1;
			for (int j = i; j < k; j++) {
				if (subMatrix[j][i] == 1) {
					pivotRow = j;
					break;
				}
			}

			if (pivotRow == -1) {
				// No 1's find the matrix is non invertible
				return false;
			}

			// swap lignes
			int[] temp = subMatrix[i];
			subMatrix[i] = subMatrix[pivotRow];
			subMatrix[pivotRow] = temp;

			// elimination
			for (int j = i + 1; j < k; j++) {
				if (subMatrix[j][i] == 1) {
					for (int l = i; l < k; l++) {
						subMatrix[j][l] ^= subMatrix[i][l];
					}
				}
			}
		}

		// check for 1's
		for (int i = 0; i < k; i++) {
			if (subMatrix[i][i] != 1) {
				return false;
			}
		}

		return true;
	}

	// here we pass from the standard format to the decimale one

	public static int[] generatorMatricesFromStandardFormat(int[][][] matrices) {
		int dim = matrices.length;
		int numRows = matrices[0].length;
		int numCols = matrices[0][0].length;

		int[] genMat = new int[dim * numCols];
		int r, c, j; // Row r, column c, dimension j.
		for (j = 0; j < dim; ++j) {
			for (r = 0; r < numRows; ++r) {
				for (c = 0; c < numCols; ++c) {
					if (matrices[j][r][c] > 0)
						genMat[j * numCols + c] += (1 << (numRows - 1 - r));
				}
			}
		}
		return genMat;
	}

}