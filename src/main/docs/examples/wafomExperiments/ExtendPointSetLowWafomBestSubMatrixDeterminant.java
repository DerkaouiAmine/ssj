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
 * This algorithm for searching low WAFOM was presented by Shin Harase. In:
 * Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357. It presents an
 * efficient method to obtain a low-WAFOM. It proceeds gradually and extends the
 * point set \(P_d = \{\mathbf{x}_0, \dots, \mathbf{x}_{2^d-1}\} \supsetneq
 * P_{d-1}\) for \(0 \leq d \leq k\), where each \(P_d\) is a subset of \(P\)
 * and \(|P|=2^k\). The algorithm tests for each \(d\) a large number of
 * different point sets, keeping the best ones in terms of WAFOM.
 * 
 * This class is similar to ExtendPointSetLowWafomBestColumn, but here our focus
 * is directly on constructing the columns of the generating matrix C bit by bit
 * at random. The objective of this class is to determine the efficiency to find
 * an invertible generating a matrix purely through randomization. Due to the
 * determinant calculations, this class operates slowly. For a faster
 * alternative, refer to ExtendPointSetLowWafomBestSubMatrixRank. The results in
 * PDF come from the class ExtendPointSetLowWafomBestSubMatrixRank
 */

public class ExtendPointSetLowWafomBestSubMatrixDeterminant {// je prend par colonne
	public static void main(String[] args) throws FileNotFoundException {
		long startTime = System.currentTimeMillis();
		int dim = 5;
		int k = 22;
		// We set w=30 to be able to use the WAFOM Accelerate methods with w=30 and q=3.
		int w = 30;

		/**
		 * This refers to the number of randomizations assigned to each column of the
		 * generator matrix.
		 */
		int nbreTest = 7000;
		double[] bestWAFOMS = new double[k];
		double[] worstWAFOMS = new double[k];
		int[][] bestColonnes = new int[k][w * dim];
		int[][] worstColonnes = new int[k][dim];

		int[] bestcolumn = new int[w * dim];

		int[] BestGeneratorMatrices = new int[dim * k];

		int[][] BEstGeneramatrix = new int[k][dim];

		/**
		 * This is the generator matrix that we will modify as we progress with our
		 * tests.
		 */

		int[][] Generamatrix = new int[dim][0];
		int[][][] optimaleMatrix = new int[dim][w][k];

		// We take point sets from 2 points to 2^23 points.
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

				// Check the invertibility of each 2D matrix in 'optimaleMatrix'.
				boolean allInvertible = true;
				for (int j = 0; j < optimaleMatrix.length; j++) {
					if (!isMatrixInvertible1(optimaleMatrix[j])) {
						allInvertible = false;
						break;
					}
				}

				// If any matrix is not invertible, continue to the next iteration.
				// Otherwise, exit the loop as an invertible matrix has been found.
				if (!allInvertible) {
					continue;

				} else {
					i1 = nbreTest;// Here we keep the first invertible matrix
				}
				int index = 0;
				for (int j = 0; j < dim; j++) {
					int[] colonne1 = extraireColonne(optimaleMatrix[j], d);
					for (int i = 0; i < colonne1.length; i++) {
						bestcolumn[i + index] = colonne1[i];

					}
					index += w;

				}

				int[] generate = generatorMatricesFromStandardFormat(optimaleMatrix);

				DigitalNetBase2 dn = new DigitalNetBase2(d + 1, w, dim, generate);

				/**
				 * We calculate the WAFOM each time
				 */
				WafomFast lowWaf = new WafomFast(dn, w, 3);

				double currentWafom = lowWaf.calcWafom();

				/**
				 * We track the best and worst WAFOM values along with their corresponding best
				 * and worst columns of matrix L, and the generating matrices.
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

			System.out.println("Best Generatrice Matrix");
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
					L[j][l][c] = stream.nextInt(0, b - 1);
				}
			}
		}
		return L;
	}

	// Here, we fix the previous "d-1" columns and generate the d-th column
	// randomly.
	public static int[][][] generateMatrix1(int s, int w, int k, int b, int[][][] r, RandomStream stream) {
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
				L[j][l][k - 1] = stream.nextInt(0, b - 1);
			}
		}

		return L;
	}

	//////////////////////////////////////////////////////////////////////////

	public static boolean isMatrixInvertible(int[][] matrix) {
		// Calculez le rang de la matrice
		int rank = calculateRank(matrix);

		// Une matrice est inversible si son rang est égal à k
		return rank == matrix[0].length;
	}

	public static boolean isMatrixInvertible1(int[][] matrix) {
		int det = calculateDeterminant(matrix);
		return det != 0;
	}

	// Gaussian Elimination to get the rank
	public static int calculateRank(int[][] matrix) {
		int rowCount = matrix.length;
		int colCount = matrix[0].length;

		int rank = Math.min(rowCount, colCount);

		// Convert the matrix to its reduced row echelon form
		for (int row = 0; row < rank; row++) {
			// If the diagonal value is 0, swap the row with another one
			if (matrix[row][row] == 0) {
				boolean reduced = false;
				for (int i = row + 1; i < rowCount; i++) {
					if (matrix[i][row] != 0) {
						swapRows(matrix, row, i);
						reduced = true;
						break;
					}
				}
				if (!reduced) {
					rank--;
					matrix[row][row] = -1; // Mark this entry as processed for subsequent steps
				}
			}

			// Make the other elements in this column 0
			for (int i = 0; i < rowCount; i++) {
				if (i != row) {
					int factor = matrix[i][row] / matrix[row][row];
					for (int j = row; j < colCount; j++) {
						matrix[i][j] -= factor * matrix[row][j];
					}
				}
			}
		}

		return rank;
	}

	/**
	 * Swaps two rows of a given matrix.
	 * 
	 * @param matrix The matrix containing the rows.
	 * @param row1   The index of the first row to swap.
	 * @param row2   The index of the second row to swap.
	 */
	public static void swapRows(int[][] matrix, int row1, int row2) {
		int[] temp = matrix[row1];
		matrix[row1] = matrix[row2];
		matrix[row2] = temp;
	}

//here we pass from the standard format to the decimale one
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

	/**
	 * Computes the determinant of a given square matrix using Laplace expansion.
	 * 
	 * @param matrix The input matrix to compute the determinant for.
	 * @return The determinant of the matrix.
	 */

	public static int calculateDeterminant(int[][] matrix) {
		int n = matrix[0].length; // Nombre de colonnes (k)

		// take only the submatrix k x k
		int[][] matrix1 = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				matrix1[i][j] = matrix[i][j];
			}

		}

		// Cas de base pour une matrice 1x1
		if (n == 1) {
			return matrix1[0][0];
		}

		// Initialze
		int det = 0;

		// Iterate over the first row of the matrix
		for (int i = 0; i < n; i++) {
			// Compute the determinant of the submatrix (by excluding the first row and the
			// ith column)
			int[][] submatrix = new int[n - 1][n - 1];
			for (int j = 1; j < n; j++) {
				for (int k = 0; k < n; k++) {
					if (k < i) {
						submatrix[j - 1][k] = matrix1[j][k];
					} else if (k > i) {
						submatrix[j - 1][k - 1] = matrix1[j][k];
					}
				}
			}

			// Calculer the submatrix determinant
			int subDet = calculateDeterminant(submatrix);

			// Adjust the main determinant based on the current element and the alternating
			// sign
			det ^= matrix[0][i] & ((i % 2 == 0) ? 1 : -1) * subDet;
		}

		return det;
	}

}