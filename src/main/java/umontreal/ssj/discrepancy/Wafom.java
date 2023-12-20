package umontreal.ssj.discrepancy;

import umontreal.ssj.hups.DigitalNetBase2;

/**
 * This class computes the Walsh Figure of Merit (WAFOM) for a Digital Net Base 2.
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
 * An alternative to this measure is introduced by Akehito Yoshiki. Bounds on Walsh coefficients
 * by dyadic difference and a new Koksma-Hlawka type inequality for Quasi-Monte Carlo integration. 2015.
 * It follows the same formula as the original WAFOM but with (j+1) as the weight. The formula is:
 *
 * WAFOM(P) = \frac{1}{|P|} \sum_{B \in P} \left\{ \prod_{i=1}^{s} \prod_{j=1}^{w} \left([1 + (-1)^{(-1)^{x_{i,j}}}2^{-(j+1)}] - 1\right)\right\}
 *
 * Another alternative specific to the shifted digital net, described in T. Goda et al.
 * “The Mean Square Quasi-Monte Carlo Error for Digitally Shifted Digital Nets”. In:
 * Monte Carlo and Quasi-Monte Carlo Methods 2014. Ed. by R. Cools and D. Nuyens.
 * Berlin: Springer-Verlag, 2016, pp. 331–350, proposes a similar criterion. By taking the square root, we have:
 *
 * W(P) = \sqrt{\frac{1}{|P|} \sum_{B \in P} \left\{ \prod_{i=1}^{s} \prod_{j=1}^{w} \left([1 + (-1)^{(-1)^{x_{i,j}}}2^{-j}] - 1\right)\right\}}
 * In this case, we have to take the square root of the results coming from this class. Otherwise, use WafomRQMC.
 *
 * We have defined a general case in this class. Setting c=1 gives the original WAFOM,
 * c=2 gives the WAFOM of Yoshiki, and factor=2 gives the criterion proposed by Goda.
 */


public class Wafom {

	/**
	 * The 'factor' here is defined as 1 for the Original definition of the WAFOM. For the definition provided by Goda, we need to set 'factor' to 2
	 * and take the square root of the final result.
	 */
	private static final double factor = 1.0;//=2 dans un article pour WAFOm equivalent RQMC
	private static DigitalNetBase2 dn;

	/**
	 * For the original definition of WAFOM, set c=1. For Yoshiki's definition, use c=2.
	 */
	private static double c=1;
	private static int s;//dimension
	private static int w;//precision
	/**
	 * Constructor for a Digital Net Base 2. The parameter 'c' is provided according to the desired method,
	 * 's' represents the dimension, 'w' is the precision, and 'N = 2^k' is the number of points.
	 *
	 * @param digitalNet The Digital Net Base 2.
	 * @param c The parameter 'c' according to the chosen method.
	 * @param s The dimension.
	 * @param w The precision.
	 * @param k The exponent determining the number of points (N = 2^k).
	 */
	public Wafom (DigitalNetBase2 dn, double c, int s,int w,int k) {
		this.dn=dn;
		this.c=c;
		this.s=s;
		this.w=w;

	}

	/**
	 * Here, we compute -1^x_ij from the WAFOM formula, or simply use the equivalent expression 1 - 2 * x_ij for a faster computation.
	 */

	private static int m1p(int x) {
		return 1 - 2 * x;
	}


	/**
	 * This method calculates the part \(\prod_{i=1}^{s} \prod_{j=1}^{w} \left([1 + (-1)^{(-1)^{x_{i,j}}}2^{-(j)}] - 1\right)\) of the formula.
	 */
	private static double calcWafomSub(int[] point) {


		double prod = 1.0;
		int startIndex=0;

		for (int i = 0; i < s; i++) {
			int index=startIndex*w;
			for (int j = 0; j < w; j++) {
				//int bij = (int) (((long)w >> (N - j - 1)) & 1);
				int bij=point[j+index];
				//System.out.println("Index:"+(j+index)+":BIJ:"+bij);

				double exponent = -factor * (j + 2 - c);
				double result = 1.0 / (1L << (int)(-exponent));
				//double p = 1.0 + m1p(bij) * Math.pow(2.0, -factor * (j + 2 - c));

				double p = 1.0 + m1p(bij) * result;

				prod *= p;
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

			double sub = calcWafomSub(point);
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
			double sub = calcWafomSub(point);
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