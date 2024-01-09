package wafomExperiments;

import java.io.*;
import java.util.Arrays;

import umontreal.ssj.hups.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.mcqmctools.*;

/**
 * In this class, we test our best point sets found with the WAFOM using the Genz functions.
 * 
 * Here, \( \mathbf{a} = (a_1, \dots, a_s) \) represents the degrees of difficulty, and \( \mathbf{u} = (u_1, \dots, u_s) \) are the shift parameters.
 * The values of \( \mathbf{a} \) and \( \mathbf{u} \) are chosen as random uniform vectors in \( [0,1]^s \), subject to the constraint \( \sum_{i=1}^s a_i = h_j \).
 * The values of \( h_j \) are set differently depending on the function.
 * They are taken from V. Barthelmann, E. Novak, and K. Ritter. “High dimensional polynomial interpolation on sparse grids”.
 * For this, we generate several (here 20) different sets for \( \mathbf{a} \) and \( \mathbf{u} \), resulting in several (here 20) examples for each function by changing \( \mathbf{a} \) and \( \mathbf{u} \).
 * For each sample \( |P| = 2^d \) where \( 1 \leq d \leq k \), we calculate the variance.
 * In the end, we take the median of the obtained variances for all the sets \( \mathbf{a} \) and \( \mathbf{u} \).
 */

/**
 * Lines sum=genz_functions.oscillatory(x, a, u); and double hi = h1[0]; need to
 * be changed for testing with other Genz functions.
 *
 */

public class GenzfunctionTestWithBestPointSets implements MonteCarloModelDouble {
	static int dim;
	double sum;
	int nbreSet;
	int w;
	static double[] a;
	int k;
	static double[] u;

	/**
	 * Constructor for GenzfunctionTestWithBestPointSets class.
	 *
	 * @param dim     Dimensionality of the function.
	 * @param w       Precision.
	 * @param k       Parameter.
	 * @param nbreSet Number of sets to generate.
	 * @param a       Array representing the degrees of difficulty.
	 * @param u       Array representing the shift parameters.
	 */
	public GenzfunctionTestWithBestPointSets(int dim, int w, int k, int nbreSet, double[] a, double[] u) {
		this.dim = dim;
		this.nbreSet = nbreSet;
		this.w = w;
		this.k = k;
		this.a = a;
		this.u = u;
	}

	public void simulate(RandomStream stream) {
		sum = 0.0;
		double[] x = new double[dim];
		for (int j = 0; j < dim; j++) {
			x[j] = stream.nextDouble();

		}

		/**
		 * We need to make changes to obtain the desired function. Also, modify the
		 * line: double hi = h1[i]; according to the desired method 'i'.
		 */
		// sum=genz_functions.gaussian(x, a, u);
		// sum=genz_functions.oscillatory(x, a, u);//changer le h
		sum = genz_functions.discontinuous(x, a, u);
	}

	public double getPerformance() {
		return sum;
	}

	public String toString() {
		return "Test function for MC and RQMC: ";
	}

	public static void main(String[] args) throws IOException {
		int k = 21;
		int w = 30;
		int dim = 5;

		Tally stats = new Tally("Simul MC");
		RandomStream random = new MRG32k3a();
		int nbreSet = 20;

		double[][] uTotal = new double[nbreSet][dim];

		/**
		 * hj values: For s=10, OSc=h0, productPeak=h1, CornerPeak=h2, Gaussian=h3,
		 * Continuous=h4, Discontinuous=h5
		 */
		double[] hGeneral = { 9, 7.25, 1.85, 7.03, 20.4, 4.3 };

		double[] h1 = new double[hGeneral.length];

		/**
		 * We adapt the values according to our dimension.
		 */
		for (int i = 0; i < hGeneral.length; i++) {
			h1[i] = dim * hGeneral[i] / 10;

		}

		/**
		 * Choose according to the desired function: OSc=0, productPeak=1, CornerPeak=2,
		 * Gaussian=3, Continuous=4, Discontinuous=5.
		 */

		double hi = h1[5];

		for (int i = 0; i < nbreSet; i++) {
			random.nextArrayOfDouble(uTotal[i], 0, dim);
		}

		double[][] aTotal = generateSetForA(nbreSet, hi, dim);// Constants a_i

		System.out.println("parametres a");
		for (double[] tab : aTotal) {
			System.out.println(Arrays.toString(tab));

		}
		System.out.println("parametres U");
		for (double[] tab1 : uTotal) {
			System.out.println(Arrays.toString(tab1));

		}

		// For RQMC+WAFOM

		System.out.println("******************Methode RQMC****************");

		int s = dim;

		// ############## S = 3 #################

		/**
		 * M0: Scrambled Sobol M1: Scrambled NX M2: Sobol M3: NX M4: Extended
		 * Construction Naive M5: Extended Sobol
		 */

		// ############## S = 12 #################

		/**
		 * M2: Scrambled NX using RQMC W(P)
		 */

		// ############## S = 16 #################

		/**
		 * M0: Scrambled Sobol M1: Scrambled NX M2: Scrambled Sobol using RQMC W(P)
		 */

		/**
		 * Parameters: - Best WAFOM matrix = 1 - Worst WAFOM matrix = 0
		 */

		/**
		 * Here, adjustments are needed based on the dimension. Referring to
		 * WafomsStorage : - For s=5, choose method < 6 - For s=12, choose method < 1 -
		 * For s=16, choose method < 3 /** Note: These are the tests I've already
		 * conducted and have stored the best matrices. For any new generating matrices,
		 * they should be added to the WafomsStorage class.
		 */

		for (int methode = 0; methode < 4; methode++) {
			System.out.println("******************Methode numero" + methode + "****************");
			int BestOrWorstWaf = 0;
			/**
			 * s=12 If method == 0, increment the BestOrWorstWaf counter.
			 */

			if (methode == 2 || methode == 3 || methode == 4 || methode == 5) {
				BestOrWorstWaf++;
			}

			while (BestOrWorstWaf < 2) {
				if (BestOrWorstWaf == 0) {
					System.out.println("******************Worst****************");
				}
				if (BestOrWorstWaf == 1) {
					System.out.println("******************Best****************");
				}

				/**
				 * We retrieve our best generator matrix.
				 */

				WafomsStorage storage = new WafomsStorage(s, methode, BestOrWorstWaf);
				int[][] Generamatrix = storage.getWafomsTabl();

				Tally stats1 = new Tally("Simul RQMC");

				int replicates = 100;

				double[] moyenneRQMC = new double[nbreSet];
				double[] varianceRQMC = new double[nbreSet];
				double[] medianeDesVariance = new double[Generamatrix.length];

				double[][] moyennesRQMCGenerale = new double[Generamatrix.length][nbreSet];
				double[][] varianceRQMCGenerale = new double[Generamatrix.length][nbreSet];

				for (int k1 = 1; k1 <= Generamatrix.length; k1++) {
					/**
					 * We create our best point set.
					 */
					DigitalNetBase2 dn = new DigitalNetBase2(k1, w, dim, Generamatrix[k1 - 1]);

					for (int i = 0; i < nbreSet; i++) {

						stats1.init();

						PointSetRandomization rand = new RandomShift(new MRG32k3a());

						RQMCExperiment.simulReplicatesRQMC(
								new GenzfunctionTestWithBestPointSets(dim, w, k, nbreSet, aTotal[i], uTotal[i]), dn,
								rand, replicates, stats1);

						moyenneRQMC[i] = stats1.average();
						varianceRQMC[i] = stats1.variance();

					}
					moyennesRQMCGenerale[k1 - 1] = moyenneRQMC;
					varianceRQMCGenerale[k1 - 1] = varianceRQMC;

					/**
					 * We retrieve the median of our variance.
					 */

					medianeDesVariance[k1 - 1] = calculateMedian(varianceRQMC);
					System.out.println(medianeDesVariance[k1 - 1]);

				}
				BestOrWorstWaf++;
			}

		}

	}

	/**
	 * Calculate the median of a vector.
	 */
	public static double calculateMedian(double[] arr) {
		Arrays.sort(arr);
		int n = arr.length;

		if (n % 2 == 0) {
			double middle1 = arr[n / 2 - 1];
			double middle2 = arr[n / 2];
			return (middle1 + middle2) / 2.0;
		} else {
			return arr[n / 2];
		}
	}

	/**
	 * We generate for our several a each one in a single ligne.
	 */

	public static double[][] generateSetForA(int numSets, double targetSum, int dim) {
		double[][] aTotal = new double[numSets][dim];
		RandomStream random = new MRG32k3a();

		for (int i = 0; i < numSets; i++) {
			double sum = 0.0;
			for (int j = 0; j < dim - 1; j++) {
				double aj = random.nextDouble();
				aTotal[i][j] = aj;
				sum += aj;
			}

			double aLast = targetSum - sum;
			if (aLast < 0 || aLast > 1) {
				/**
				 * If the constraint is not met, redo the generation.
				 */
				i--;
				continue;
			}

			aTotal[i][dim - 1] = aLast;
		}
		return aTotal;
	}

	/**
	 * This method is essential for replacing our generator matrix with a vector and
	 * using this vector to create our Sobol point.
	 */

	public static int[] matrixToVector(int[][] matrix) {
		if (matrix == null || matrix.length == 0 || matrix[0].length == 0) {
			return new int[0];
		}

		int m = matrix.length;
		int k = matrix[0].length;
		int[] vector = new int[m * k];

		int index = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < k; j++) {
				vector[index++] = matrix[i][j];
			}
		}

		return vector;
	}

	/**
	 * This method is essential for constructing the generator matrix because each
	 * time we extract either a single element or the preceding submatrix, which
	 * will already be optimal.
	 */

	public static int[][] extractSubMatrix(int[][] originalMatrix, int endRow, int endCol) {
		int startCol = 0;
		int startRow = 0;
		// Vérifier les limites des indices pour éviter les erreurs
		int numRows = originalMatrix.length;
		int numCols = originalMatrix[0].length;

		if (startRow < 0 || endRow >= numRows || startCol < 0 || endCol >= numCols) {
			throw new IllegalArgumentException("Les indices de sous-matrice sont hors limites.");
		}

		// Calculer la taille de la sous-matrice
		int subMatrixRows = endRow - startRow + 1;
		int subMatrixCols = endCol - startCol + 1;

		// Créer la sous-matrice
		int[][] subMatrix = new int[subMatrixRows][subMatrixCols];

		// Copier les éléments de la matrice originale dans la sous-matrice
		for (int i = startRow; i <= endRow; i++) {
			for (int j = startCol; j <= endCol; j++) {
				subMatrix[i - startRow][j - startCol] = originalMatrix[i][j];
			}
		}

		return subMatrix;
	}

}
