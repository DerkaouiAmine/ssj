package umontreal.ssj.discrepancy;

import java.io.FileNotFoundException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.stat.Tally;



/**
 * This class computes another figure of merit close to the Walsh Figure of Merit (WAFOM)
 * but specific to the shifted digital net, as described in T. Goda et al.
 * “The Mean Square Quasi-Monte Carlo Error for Digitally Shifted Digital Nets”. In:
 * Monte Carlo and Quasi-Monte Carlo Methods 2014. Ed. by R. Cools and D. Nuyens.
 * Berlin: Springer-Verlag, 2016, pp. 331–350. It has the following form:
 *
 * W(P) = \sqrt{\frac{1}{|P|} \sum_{B \in P} \left\{ \prod_{i=1}^{s} \prod_{j=1}^{w}
 * \left([1 + (-1)^{(-1)^{x_{i,j}}}2^{-2*j}] - 1\right)\right\}}
 *
 * where |P| is the number of points in P, s is the dimension of the points,
 * w is the precision, x_{i,j} is the j-th coordinate of point x in the i-th dimension.
 * 
 *  An alternative to this measure is introduced by Akehito Yoshiki. Bounds on Walsh coefficients
 * by dyadic difference and a new Koksma-Hlawka type inequality for Quasi-Monte Carlo integration. 2015.
 * It follows the same formula as the original WAFOM but with (j+1) as the weight. It has the following form:
 *
 * W(P) = \sqrt{\frac{1}{|P|} \sum_{B \in P} \left\{ \prod_{i=1}^{s} \prod_{j=1}^{w}
 * \left([1 + (-1)^{(-1)^{x_{i,j}}}2^{-2*(j+1)}] - 1\right)\right\}}
 *
 * This discrepancy measure, specific to digitally shifted nets, is designed and shows
 * excellent results for alpha-smooth functions (see J. Dick. “On Quasi-Monte Carlo Rules
 * Achieving Higher Order Convergence”. In: Monte Carlo and Quasi-Monte Carlo Methods 2008).
 */


public class WafomRQMC {

	/**
	 * We kept the definition the same as for WAFOM, but here, factor=2 for this method, compared to 1 for WAFOM.
	 */
	private static final double factor = 2.0;
	private static DigitalNetBase2 dn;
	private static double c=1;//by default and 2 if we want the Yoshiki definition
	private static int s;
	private static int w;
	private static int k;
	private static PointSetRandomization rand;

	/**
	 * Constructor for a Digital Net Base 2. The parameter 'c' is provided according to the desired method,
	 * 's' represents the dimension, 'w' is the precision, and 'N = 2^k' is the number of points.
	 *
	 * @param digitalNet The Digital Net Base 2.
	 * @param c The parameter 'c' according to the chosen method.
	 * @param s The dimension.
	 * @param w The precision.
	 * @param k The exponent determining the number of points (N = 2^k).
	 * @param rand is essential to establish a digital shift	
	 * 
	 */

	public WafomRQMC (DigitalNetBase2 dn, double c, int s,int w,int k,  RandomShift rand) {
		this.dn=dn;
		this.c=c;
		this.s=s;
		this.w=w;
		this.k=k;
		this.rand=rand;
	}
	/**
	 * Here, we compute -1^x_ij from the WAFOM formula, or simply use the equivalent expression 1 - 2 * x_ij for a faster computation.
	 */
	private static int m1p(int x) {
		return 1 - 2 * x;
	}


	/**
	 * This method calculates the part \(\prod_{i=1}^{s} \prod_{j=1}^{w} \left([1 + (-1)^{(-1)^{x_{i,j}}}2^{-2*(j)}] - 1\right)\) of the formula.
	 */
	private static double calcWafomSub(int[] point) {
		int N=w;

		double prod = 1.0;
		int startIndex=0;

		for (int i = 0; i < s; i++) {

			int index=startIndex*w;
			for (int j = 0; j < w; j++) {

				int bij=point[j+index];

				double p = 1.0 + m1p(bij) * Math.pow(2.0, -factor*(j+2-c));
				prod *= p;
			}
			startIndex++;
		}

		return prod ;
	}



	/**
	 * In this method, we apply the previous method calcWafomSub for each point among the N = 2^k points.
	 * Simply sum the results and divide by |P| to obtain the WAFOM for our point set.
	 */
	public static double calcWafom(){
		double sum = 0.0;

		rand.randomize(dn);

		// dn.addRandomShift();

		double[][] tab= dn.formatPointsTab();

		int [][]pointSet=convertDecimalToBinary(tab,w);
		long num = pointSet.length;
		for (long i = 0; i < num; i++) {
			int[] point = pointSet[(int) i];

			double sub = calcWafomSub(point);

			sum += sub;
		}

		return Math.sqrt(-1+(sum/num));
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