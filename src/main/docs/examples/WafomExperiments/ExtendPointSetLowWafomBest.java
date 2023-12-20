package WafomExperiments;



import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;

/**
 * This algorithm for searching low WAFOM was presented by Shin Harase.
 * In: Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357.
 * It presents an efficient method to obtain a low-WAFOM. It proceeds gradually
 * and extends the point set \(P_d = \{\mathbf{x}_0, \dots, \mathbf{x}_{2^d-1}\} \supsetneq P_{d-1}\) for \(0 \leq d \leq k\),
 * where each \(P_d\) is a subset of \(P\) and \(|P|=2^k\). The algorithm tests for each \(d\) a large number of different point sets,
 * keeping the best ones in terms of WAFOM.
 * 
 * We will work on the generator matrices, so all randomizations are done on the columns of the generator matrices gradually for each d.
 *  */

public class ExtendPointSetLowWafomBest {//je prend par colonne
	public static void main(String[] args) throws FileNotFoundException { 
		long startTime = System.currentTimeMillis();

		int dim=16;
		int k=23;


		// We set w=30 to be able to use the WAFOM Accelerate methods with w=30 and q=3.
		int w=30;


		/**
		 * This refers to the number of randomizations assigned to each column of the generator matrix.
		 */
		int nbreTest=3000;


		double[]bestWAFOMS=new double[k];
		double[]worstWAFOMS=new double[k];
		int[][] bestColonnes=new int [ k][dim];
		int[][] worstColonnes=new int [ k][dim];




		/**
		 * This is the generator matrix that we will modify as we progress with our tests.
		 */


		int[][] Generamatrix = new int[dim][0]; 

		// We take point sets from 2 points to 2^23 points.

		for(int d=0;d<k;d++) {


			System.out.println("Je calcule pour k="+(d+1));

			double smallestWafom=(1<<k);
			double biggestWafom=-1;


			// For each k, we test 'nbreTest' times.

			for(int i1=0;i1<nbreTest;i1++) {
				DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
				Sob.leftMatrixScramble(new MRG32k3a());

				int [][]	tableauPrincipalGenerLMS =  Sob.getGeneratorMatrices(k);//ici sa sera un i au lieu de 0





				int[] colonne1 = extraireColonne(tableauPrincipalGenerLMS, d);//ici i1 a la place de 0

				if (i1==0) {
					Generamatrix = concatenateMatrixAndColumn(Generamatrix, colonne1);}
				else {
					Generamatrix=replaceColumn(Generamatrix,colonne1, d);
				}


				int vect[]=matrixToVector(Generamatrix);


				/**
				 * For each "d," we generate columns of the generator matrix 'nbreTest' times.
				 */
				DigitalNetBase2 Sob1 = new SobolSequence(d+1,w,dim,vect);

				double [][]SobPoints= Sob1.formatPointsTab();








				/**
				 * We calculate the WAFOM each time 
				 */


				Wafom lowWaf=new Wafom (Sob1,1,dim,w,d+1);

				double currentWafom =lowWaf.calcWafom1(SobPoints,d+1);

				/**
				 * We keep track of the best and worst WAFOM values.
				 */

				if (currentWafom < smallestWafom) {
					smallestWafom = currentWafom;

					bestWAFOMS[d]=smallestWafom;
					bestColonnes[d]=colonne1;

				}

				if (currentWafom > biggestWafom) {
					biggestWafom= currentWafom;

					worstWAFOMS[d]=biggestWafom;
					worstColonnes[d]=colonne1;

				}













			}
			long endTime = System.currentTimeMillis();
			long executionTime = endTime - startTime;

			System.out.println(executionTime);

			/**
			 * We keep the generator matrix with the best column that generates the best points in terms of WAFOM
			 * to extend our set to a larger d.
			 */
			Generamatrix=replaceColumn(Generamatrix,bestColonnes[d], d);

			System.out.println("best WAFOMS");
			System.out.println(Arrays.toString( bestWAFOMS));


			System.out.println("BEST COLOMNS");
			for (int i = 0; i <  bestColonnes[0].length; i++) {
				for (int j = 0; j <  bestColonnes.length; j++) {
					System.out.print(bestColonnes[j][i] + " ");
				}
				System.out.println();
			}


			System.out.println("worst WAFOMS");
			System.out.println(Arrays.toString( worstWAFOMS));


			System.out.println("worst COLOMNS");
			for (int i = 0; i <  worstColonnes[0].length; i++) {
				for (int j = 0; j <  worstColonnes.length; j++) {
					System.out.print(worstColonnes[j][i] + " ");
				}
				System.out.println();
			}

		}



	}



	/**
	 * This method is essential for constructing our general generator matrix, starting with an empty matrix.
	 * We call this method each time to add the best column in terms of WAFOM to this matrix.
	 */

	public static int[][] concatenateMatrixAndColumn(int[][] matrix, int[] colonne1) {
		int numRows = matrix.length;
		int numCols = matrix[0].length;

		if (numRows != colonne1.length) {
			throw new IllegalArgumentException(" The number of rows in the matrix must be equal to the size of the column vector.");
		}

		int[][] resultMatrix = new int[numRows][numCols + 1];

		for (int i = 0; i < numRows; i++) {
			for (int j = 0; j < numCols; j++) {
				resultMatrix[i][j] = matrix[i][j];
			}
			resultMatrix[i][numCols] = colonne1[i];
		}

		return resultMatrix;
	}


	/**
	 * This method is essential because we generate generator matrices with k columns,
	 * and we extract just the k-th column 'nbreTest' times to test it.
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
	 * When we extract our 'nbre sets' of columns, we should replace them in our generator matrix
	 * to test them in terms of WAFOM.
	 */

	public static int[][] replaceColumn(int[][] matrix, int[] newColumn, int columnIndex) {
		if (columnIndex < 0 || columnIndex >= matrix[0].length) {
			throw new IllegalArgumentException("L'index de colonne spécifié est invalide.");
		}

		if (newColumn.length != matrix.length) {
			throw new IllegalArgumentException("Le nouveau vecteur colonne doit avoir la même taille que le nombre de lignes de la matrice.");
		}

		for (int i = 0; i < matrix.length; i++) {
			matrix[i][columnIndex] = newColumn[i];
		}

		return matrix;
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



	/**
	 * This does the reverse of the previous method if needed.
	 */
	public static int[][] vectorToMatrix(int[] vector) {
		int numRows = vector.length;
		int[][] matrix = new int[numRows][1];

		for (int i = 0; i < numRows; i++) {
			matrix[i][0] = vector[i];
		}

		return matrix;
	}





}