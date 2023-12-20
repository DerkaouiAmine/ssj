package umontreal.ssj.discrepancy;

import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;

/**
 * Here, we adopt the same construction definition as proposed for the class WAFOM.
 * 
 *  This class computes the Walsh Figure of Merit (WAFOM) for a Digital Net Base 2.
 * @cite {M. Matsumoto, M. Saito, and K. Matoba. “A Computable Figure of Merit for
 * Quasi-Monte Carlo Point Sets”. In: Mathematics and Computers in Simulation
 * 83.287 (2014), pp. 1233–1250}
 * It is given by:
 *
 * WAFOM(P) = \frac{1}{|P|} \sum_{B \in P} \left\{ \prod_{i=1}^{s} \prod_{j=1}^{w} \left([1 + (-1)^{(-1)^{x_{i,j}}}2^{-j}] - 1\right)\right\}
 *
 * where |P| is the number of points in P, s is the dimension of the points,
 * w is the precision, x_{i,j} is the j-th coordinate of point x in dimension i-th dimension.
 *
 * This discrepancy measure specific to digital nets is designed and shows excellent results
 * for alpha-smooth functions (see J. Dick. “On Quasi-Monte Carlo Rules Achieving Higher Order Convergence”.
 * In: Monte Carlo and Quasi-Monte Carlo Methods 2008).
 * 
 * 
 * 
 * The difference here lies in the method proposed by Shin Harase in "A search for extensible low-WAFOM point sets".
 * In: Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357. The approach involves taking a coordinate of a point
 * with its 'w' precision and dividing it into 'q' subparts, each containing 'l' elements. Thus, a coordinate of a point
 * is expressed as X^{i} = d_1^i, ..., d_q^i, where for 1 <= c <= q, our d_c^i is x_{i, (c-1)*l-1}, ..., x_{i, c*l}.
 * Instead of computing \prod_{j=1}^{w} \left([1 + (-1)^{(-1)^{x_{i,j}}}2^{-j}] - 1, we replace it with
 * \prod_{1<=c<=q} table_c[d_c^i], where the table_c[d_c^i] are precalculated and take the form
 * table_c[d_c^i] = \prod_{1<=j<=l}(1+(-1)^{x_{i,(c-1)*l+j}} * 2^{-((c-1)*l+j+1)}). This significantly reduces computation time.
 */


public class WafomAcceler {
	private static DigitalNetBase2 dn;
	private static int dim;//dimension s
	private static int w;
	private  static int q;

	/**
	 * Constructor for a Digital Net Base 2. The parameter 'c' is provided according to the desired method,
	 * 's' represents the dimension, 'w' is the precision, and 'N = 2^k' is the number of points.
	 *
	 * @param digitalNet The Digital Net Base 2.
	 * @param c The parameter 'c' according to the chosen method.
	 * @param s The dimension.
	 * @param w The precision.
	 * @param k The exponent determining the number of points (N = 2^k).
	 * @param q Instead of taking all the 'w' bits, we divide it into 'q' equal parts.

	 * 
	 */

	public WafomAcceler (DigitalNetBase2 dn, double c, int dim,int w,int k,int q) {
		this.dn=dn;
		this.dim=dim;
		this.w=w;
		this.q=q;



	}


	/**
	 * Here, we compute -1^x_ij from the WAFOM formula, or simply use the equivalent expression 1 - 2 * x_ij for a faster computation.
	 */
	private static int m1p(int x) {
		return 1 - 2 * x;
	}





	/**
	 * Here, for each point in the 's' dimensions, we precalculate its table_c[d_c^i].
	 */

	private static double[] calctabl(int[] point) {

		/**
		 * Generates an error because we need to divide 'w' into 'q' equal parts, and therefore, the division must result in an integer.
		 */
		if (w % q != 0)   {throw new ArithmeticException("The division of w by q is not an integer.");}


		/**
		 * If q=0, an error is returned.
		 */
		if (q == 0)  { System.err.println("Warning: q is equal to zero, division by zero may occur.");}

		int l=w/q;

		double[]tabledci=new double [q*dim];

		double prod=1;
		int startIndex=0;
		int index1=0;

		for (int i = 0; i < dim; i++) {
			//double w1 = point[i];
			int index=startIndex*w;

			for (int c = 1; c <=q; c++) {
				prod = 1.0;

				for (int j = 1; j <= l; j++)
				{

					int bij=point[(c-1)*l+(j-1)+index];
					double exponent = -((c-1)*l+(j-1)+1);
					double result = 1.0 / (1L << (int)(-exponent));
					/*
					 * we do not use result Math.pow(2.0, -((c-1)*l+(j-1)+1)) for a faster computation
					 * */
					prod*=1.0 + m1p(bij) *result;


				}

				tabledci[index1]=prod;
				index1++;

			} 
			startIndex++;

		}


		return tabledci;

	}




	/**
	 * This method calculates the part \prod_{1<=c<=q} table_c[d_c^i] of the formula.
	 */
	private static double calcWafomSub(double [] dci) {
		double prod = 1.0;
		int startIndex=0;
		int indexxx=0;

		for (int i = 0; i < dim; i++) {
			for (int c = 1; c <=q ; c++) {
				//prod*=dci[(c-1)+index];

				double p =dci[indexxx];

				prod*=p;
				indexxx++;
			}

			startIndex++;
		}

		return prod - 1.0;
	}




	/**
	 * In this method, we apply the previous method calcWafomSub for each point among the N = 2^k points.
	 * Simply sum the results and divide by |P| to obtain the WAFOM for our point set.
	 */



	public static double calcWafom(){





		double sum = 0.0;

		double[][] tab= dn.formatPointsTab();

		int [][]pointSet=convertDecimalToBinary(tab,w);
		long num = pointSet.length;
		for (long i = 0; i < num; i++) {
			int[] point = pointSet[(int) i];
			double []tabledc=calctabl( point);

			double sub = calcWafomSub(tabledc);
			sum += sub;
		}

		return sum/ num;
	}


	/**
	 * Here is the same definition of WAFOM, but taking an array instead of a digital net.
	 * I chose not to use this as a constructor for a more general usage.
	 */

	public static double calcWafom1(double[][] array,int k1) {
		double sum = 0.0;
		int num1=(1<<k1);
		int[][] pointSet = convertDecimalToBinary(array, w);
		long num = pointSet.length;
		for (long i = 0; i < num; i++) {

			int[] point = pointSet[(int) i];
			double []tabledc=calctabl( point);

			double sub = calcWafomSub(tabledc);
			sum += sub;

		}
		return sum / num1;
	}



















	/**
	 * This method converts our decimal point set into a binary point set for later use in the WAFOM calculation.
	 * It takes the parameter 'decimalNumbers,' which represents our point set with 'N = 2^k' rows and 's' columns in decimal.
	 * It transforms it into base 2 with a conversion precision of 'decimalPlaces.'
	 */

	public static int[][] convertDecimalToBinary(double[][] decimalNumbers, int decimalPlaces) {
		int[][] binaryArray = new int[decimalNumbers.length][decimalNumbers[0].length * decimalPlaces];

		for (int i = 0; i < decimalNumbers.length; i++) {
			for (int j = 0; j < decimalNumbers[i].length; j++) {
				double decimalNumber = decimalNumbers[i][j];
				int[] binaryRow = new int[decimalPlaces];

				for (int k = 0; k < decimalPlaces; k++) {
					decimalNumber *= 2;
					binaryRow[k] = (int) decimalNumber;
					decimalNumber -= binaryRow[k];
				}

				System.arraycopy(binaryRow, 0, binaryArray[i], j * decimalPlaces, decimalPlaces);
			}
		}

		return binaryArray;
	}



}