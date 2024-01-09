package umontreal.ssj.discrepancy;

import java.io.FileNotFoundException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetIterator;
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
	private static final double c=1;//by default and 0 if we want the Yoshiki definition
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
	/**
	 * This construction is intentionally kept for testing purposes. 
	 * It specifically caters to cases where the number of points is less than or equal to 'k' 
	 * and when construction is done dimension by dimension.
	 * 
	 * The standard 'WafomRQMCWafomRQMC (DigitalNetBase2 dn,  RandomShift rand) as the original Wafom constructor, 
	 * can be replaced here or another constructor can be used depending on the context.
	 */

	public WafomRQMC (DigitalNetBase2 dn, int k,int w,int s,  RandomShift rand) {
		this.dn=dn;
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


		double prod = 1.0;
		int startIndex=0;

		for (int i = 0; i < s; i++) {
			int index=startIndex*w;
			for (int j = 0; j < w; j++) {
				int bij=point[j+index];
				int exponent = (int) (factor * (j + 2 - c));
				double result = 1.0 / (1L << exponent);
				double p = 1.0 + m1p(bij) * result;
				prod *= p;
			}
			startIndex++;
		}

		return (prod );
	}
	/**
	 * In this method, we apply the previous method calcWafomSub for each point among the N = 2^k points.
	 * Simply sum the results and divide by |P| to obtain the WAFOM for our point set.
	 */
	public static double calcWafom(){

		rand.randomize(dn);
		double sum = 0.0;
		long num = 1<<k;
		double[] point=new double[s];
		int [] bPoint= new int[s*w];
		int i=0;
		int expon=1<<(w+1);
		double shift=1/expon;

		PointSetIterator iter=dn.iteratorNoGray();
		PointSetIterator iterator = dn.iteratorNoGray();
		while(iter.hasNextPoint()) {
			i=0;
			int index=0;
			while(iter.hasNextCoordinate()) {
				double coorddi=iter.nextCoordinate();
				int[] currentPoint = iterator.getCachedCurPoint();//+shift
				for(int j=0;j<currentPoint.length;j++) {
					bPoint[j+index]=(int) (currentPoint[j]+shift);
				}
				index+=w;
			}	
			double sub = calcWafomSub(bPoint);
			sum += sub;
			iter.resetToNextPoint();

		}

		return Math.sqrt(-1+(sum/num));
	}







}