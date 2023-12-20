package WafomExperiments;



import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG31k3p;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

/**
 * This algorithm for searching low WAFOM was presented in my work.
 * The construction is done by constructing the general generating matrix element by element,
 * starting from the best submatrices gradually to eventually reach the complete generating matrix.
 * It presents an efficient method to obtain a low-WAFOM.
 * We will take, for each row, its decimal representation multiplied by 'w*k',
 * so we get a '1*k' vector for each dimension, resulting in 'dim*k'.
 * - We generate the first vectors \(\overline{c}_{1}^{(1)}, \dots, \overline{c}_{k}^{(1)}\) at random times.
 *   We select the vector with the smallest WAFOM each time for each \(1\leq d\leq k\).
 *
 * - For dimension \(s\) and the \(i\)-th element of the matrix, we suppose that all the column vectors
 *   \(\overline{c}^{2},\dots,\overline{c}^{s-1}\) for the previous \(s-1\) dimensions (for every \(1\leq d \leq k\))
 *   have been determined. Then we proceed element by element, fixing the \(i-1\) elements of dimension \(s\),
 *   and we generate randomly at several times the \(i\)-th element of the \(s\)-th dimension conditionally to
 *   the optimal submatrices \(s-1 \times (i)\) and all the optimal \(i-1\) elements of the \(s\)-th dimension.
 */

public class ExtendPointSetLowWafomGeneMatrixBest {
	public static void calc(int k) { 
		long startTime = System.currentTimeMillis();

		int dim=3;

		/**
		 * This refers to the number of random2izations assigned to each column of the generator matrix.
		 */
		int nbreTest=100;


		// usually we set w=30 to be able to use the WAFOM Accelerate methods with w=30 and q=3.

		int w=30;

		/**
		 * Store our best components.
		 */
		int  bestColonne=0;
		RandomStream randomStream = new MRG32k3a();





		/**
		 * This is the generator matrix that we will modify as we progress with our tests.
		 */


		int[][] Generamatrix= new int[dim][k];

		/**
		 * We will proceed dimension by dimension.
		 */
		for(int dimension=1;dimension<=dim;dimension++) {
			// We take point sets from 2 points to 2^23 points.

			for(int d=0;d<k;d++) {
				double smallestWafom=(1<<k);
				double biggestWafom=-1;
				// For each d, we test 'nbreTest' times.

				for(int i1=0;i1<nbreTest;i1++) {

					DigitalNetBase2 Sob = new SobolSequence(d+1,w,dimension);
					Sob.leftMatrixScramble(new MRG32k3a());

					int [][] generatriOriginale=Sob.getGeneratorMatrices(d+1);


					int element=generatriOriginale[dimension-1][d];


					Generamatrix[dimension-1][d]=element;

					int [][] generat=extractSubMatrix(Generamatrix,dimension-1,d);


					int vect[]=matrixToVector( generat);

					/**
					 * For each "d," we generate columns of the generator matrix 'nbreTest' times.
					 */
					 DigitalNetBase2 Sob1 = new SobolSequence(d+1,w,dimension,vect);

					 double [][]SobPoints= Sob1.formatPointsTab();

					 /**
					  * We calculate the WAFOM each time 
					  */

					 Wafom lowWaf=new Wafom (Sob1,1,dimension,w,d+1);

					 double currentWafom =lowWaf.calcWafom1(SobPoints,d+1);


					 /**
					  * We keep track of the best and worst WAFOM values.
					  */

					 if (currentWafom < smallestWafom) {
						 smallestWafom = currentWafom;


						 bestColonne=element;

					 }

					 if (currentWafom > biggestWafom) {
						 biggestWafom= currentWafom;

					 }









				}
				/**
				 * We keep the generator matrix with the best column that generates the best points in terms of WAFOM
				 * to extend our set to a larger d.
				 */
				Generamatrix[dimension-1][d]=bestColonne;


			}



		}
		long endTime = System.currentTimeMillis();
		long executionTime = endTime - startTime;

		System.out.println("time:"+executionTime);
		System.out.println("pour k="+k);
		for(int []vecteur:Generamatrix)
		{System.out.println(Arrays.toString(vecteur));}


		/**
		 * Here, due to the construction of the method, we retrieve our best WAFOM by recreating our point set
		 * and calculating the WAFOM. Otherwise, the results may be incorrect. We have observed that even if we have
		 * a very good WAFOM for one dimension, randomly adding a complete dimension with all components at once,
		 * rather than one by one, can yield poor results.
		 */


		int []vect=matrixToVector(Generamatrix);

		DigitalNetBase2 Sob = new SobolSequence(k,w,dim,vect);


		Wafom lowWaf=new Wafom (Sob,1,dim,w,k);

		double currentWafom =lowWaf.calcWafom();

		System.out.println(currentWafom);

	}

	public static void main(String[] args) throws FileNotFoundException {

		for(int k=1;k<=23;k++)
		{calc(k);}
	}

	/**
	 * This method is essential for constructing the generator matrix because each time we extract either a single element
	 * or the preceding submatrix, which will already be optimal.
	 */

	public static int[][] extractSubMatrix(int[][] originalMatrix, int endRow, int endCol) {
		int startCol = 0;
		int startRow = 0;
		// Check index boundaries to avoid errors
		int numRows = originalMatrix.length;
		int numCols = originalMatrix[0].length;

		if (startRow < 0 || endRow >= numRows || startCol < 0 || endCol >= numCols) {
			throw new IllegalArgumentException("Sub-matrix indices are out of bounds.");
		}

		// Calculate the size of the sub-matrix
		int subMatrixRows = endRow - startRow + 1;
		int subMatrixCols = endCol - startCol + 1;

		// Create the sub-matrix
		int[][] subMatrix = new int[subMatrixRows][subMatrixCols];

		// Copy elements from the original matrix into the sub-matrix
		for (int i = startRow; i <= endRow; i++) {
			for (int j = startCol; j <= endCol; j++) {
				subMatrix[i - startRow][j - startCol] = originalMatrix[i][j];
			}
		}

		return subMatrix;
	}








	/**
	 * This method is essential for replacing our generator matrix with a vector
	 * and using this vector to create our Sobol point.
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

}