package WafomExperiments;

import java.io.*;
import java.util.Arrays;

import umontreal.ssj.hups.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.mcqmctools.*;
/**
 * This class is exactly similar to GenzfunctionTestWithBestPointSets. The only difference here is that the test parameters, 'a' and 'u', are generated in advance of the experiment, and we take the median of each column to test our function.
 * Review GenzfunctionTestWithBestPointSets for more details. 
 */


/**
 * Line sum=genz_functions.oscillatory(x, a, u); and    double hi = h1[0]; need to be changed for testing with other Genz functions.
 */



public class GenzfunctionTestWithBestPointSets1 implements MonteCarloModelDouble {
	static int dim;
	double sum;
	int nbreSet;
	int w;
	static double []a;
	int k;
	static double []u;

	/**
	 * Constructor for GenzfunctionTestWithBestPointSets class.
	 *
	 * @param dim      Dimensionality of the function.
	 * @param w        Precision.
	 * @param k        Parameter.
	 * @param nbreSet  Number of sets to generate.
	 * @param a        Array representing the degrees of difficulty.
	 * @param u        Array representing the shift parameters.
	 */

	public GenzfunctionTestWithBestPointSets1  (int dim,int w,int k, int nbreSet,double[] a, double[] u) {
		this.dim = dim; 	this.nbreSet= nbreSet;
		this.w=w;			this.k=k;
		this.a=a;			this.u=u;
	}

	public void simulate (RandomStream stream) {
		sum= 0.0;
		double [] x=new double [dim];
		for (int j = 0; j < dim; j++) {
			x[j]= stream.nextDouble();

		}
		/**
		 * We need to make changes to obtain the desired function. Also, modify the line: 
		 * double hi = h1[i]; 
		 * according to the desired method 'i'.
		 */

		//sum=genz_functions.gaussian(x, a, u);
		sum=genz_functions.oscillatory(x, a, u);//changer le h
	}







	public double getPerformance () {
		return sum;
	}

	public String toString () {
		//return "Test function for MC and RQMC: Gaussian";
		return "Test function for MC and RQMC: Oscil";

	}

	public static void main(String[] args) throws IOException {
		int k=21;
		int w=30;
		int dim=3;

		Tally stats=new Tally("Simul MC");
		RandomStream random = new MRG32k3a();
		int nbreSet = 20;



		double[][] uTotal = new double[nbreSet][dim]; // Constants u_i

		/**
		 * hj values:
		 * For s=10, OSc=h0, productPeak=h1, CornerPeak=h2, Gaussian=h3, Continuous=h4, Discontinuous=h5
		 */	        double [] hGeneral= {9,7.25,1.85,7.03,20.4,4.3};


		 double [] h1= new double[hGeneral.length];//3/10*ce que a donne pour s=10


		 /**
		  * We adapt the values according to our dimension.
		  */
		 for (int i=0;i<hGeneral.length;i++) {
			 h1[i]=dim*hGeneral[i]/10;

		 }


		 /**
		  * Choose according to the desired function: OSc=0, productPeak=1, CornerPeak=2, Gaussian=3, Continuous=4, Discontinuous=5.
		  */
		 double hi = h1[0];

		 for (int i=0;i<nbreSet;i++) {
			 random.nextArrayOfDouble(uTotal[i], 0, dim);
		 }

		 double [][] aTotal= generateSetForA (nbreSet, hi, dim);// Constants a_i

		 System.out.println("les a");
		 for(double [] tab : aTotal) {
			 System.out.println(Arrays.toString(tab));

		 }
		 System.out.println("les U");
		 for(double [] tab1 : uTotal) {
			 System.out.println(Arrays.toString(tab1));

		 }


		 double[] mediansA = new double[aTotal[0].length]; // Tableau pour stocker les médianes

		 double[] mediansU = new double[uTotal[0].length]; // Tableau pour stocker les médianes

		 for (int i = 0; i < aTotal[0].length; i++) {
			 double[] column = new double[aTotal.length];
			 double [] column1=new double [uTotal.length];
			 for (int j = 0; j < aTotal.length; j++) {
				 column[j] = aTotal[j][i];
				 column1[j]=uTotal[j][i];
			 }
			 /**
			  * The difference here is that we calculate the medians before our test. Therefore, we will have only one 'a' and 'u' set with which we will test our functions.
			  */

			 mediansA[i] = calculateMedian(column);
			 mediansU[i] = calculateMedian(column1);

		 }


		 int s = dim;

		 for (int methode=0;methode<7;methode++) {
			 System.out.println("******************Methode MCVSRQMC numero"+methode+"****************");
			 int BestOrWorstWaf=0;
			 if (methode==3 || methode ==4 || methode ==5) {
				 BestOrWorstWaf++;
			 }

			 while(BestOrWorstWaf<2) {
				 if (BestOrWorstWaf==0) {
					 System.out.println("******************Worst****************");
				 }
				 if (BestOrWorstWaf==1) {
					 System.out.println("******************Best****************");
				 }
				 WafomsStorage storage = new WafomsStorage(s,methode,BestOrWorstWaf);
				 int[][] Generamatrix = storage.getWafomsTabl();


				 Tally stats1=new Tally("Simul RQMC");

				 new Tally();

				 //a quoi sa sert replicates???
				 int replicates=100;

				 for (int k1=14;k1<=Generamatrix[0].length;k1++) {
					 // System.out.println("--------------------------pour k= "+k1);




					 int [][] generat=extractSubMatrix(Generamatrix,Generamatrix.length-1,k1-1);


					 int []vect=matrixToVector(generat);

					 /**
					  * We create our best point set.
					  */


					 DigitalNetBase2 Sob = new SobolSequence(k1,w,dim,vect);
					 /**
					  * In a comparative section of my work, I used the digital shifted point set.
					  */

					 //Sob.addRandomShift(new MRG32k3a());

					 stats1.init();
					 stats.init();
					 PointSetRandomization rand = new RandomShift(new MRG32k3a());





					 System.out.println(  RQMCExperiment. makeComparisonExperimentMCvsRQMC (new GenzfunctionTestWithBestPointSets1  (dim,w, k,  nbreSet, mediansA, mediansU), new MRG32k3a(), 
							 Sob, rand, 1000000, replicates)); 





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
			double middle1 = arr[n/2 - 1];
			double middle2 = arr[n/2];
			return (middle1 + middle2) / 2.0;
		} else {
			return arr[n/2];
		}
	}



	/**
	 * We generate for  our several a each one in a single ligne.
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
	 * This method is essential for constructing the generator matrix because each time we extract either a single element
	 * or the preceding submatrix, which will already be optimal.
	 */


	public static int[][] extractSubMatrix(int[][] originalMatrix,  int endRow , int endCol) {
		int startCol=0;
		int startRow=0;
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
