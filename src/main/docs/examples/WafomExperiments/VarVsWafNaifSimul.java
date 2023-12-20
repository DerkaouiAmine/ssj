package WafomExperiments;

import java.awt.datatransfer.SystemFlavorMap;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;


/**
 * In this class, we test the relationship where, upon comparison, the points with the best WAFOM values closely correspond to those with the lowest variances.
 */




public class VarVsWafNaifSimul implements MonteCarloModelDouble {
	static int dim;
	double sum;
	int nbreSet;
	int w;
	static double []a;
	int k;
	static double []u;

	/**
	 * Constructor for the VarVsWafNaifSimul class.
	 *
	 * @param dim Dimension of the point sets.
	 * @param w Precision.
	 * @param k Number of point sets (|P| = 2^k).
	 * @param nbreSet Number of times to simulate the generation of point sets.
	 * @param a Array representing degrees of difficulty.
	 * @param u Array representing shift parameters.
	 */


	/**
	 * For more details, refer to the class GenzfunctionTestWithBestPointSets.
	 */


	public VarVsWafNaifSimul (int dim,int w,int k, int nbreSet,double[] a, double[] u) {
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

		//sum=genz_functions.gaussian(x, a, u);
		sum=genz_functions.oscillatory(x, a, u);
	}







	public double getPerformance () {
		return sum;
	}

	public String toString () {
		return "Test function for MC and RQMC: ";
	}


	public static void main(String[] args) throws FileNotFoundException { 



		long startTime = System.currentTimeMillis();
		int k=20;
		int w=30;
		int dim=3;

		Tally stats1=new Tally("Simul");
		RandomStream random = new MRG32k3a();
		int nbreSet = 20;



		double[][] uTotal = new double[nbreSet][dim]; // Constants u_i

		/**
		 * For more details, refer to the class GenzfunctionTestWithBestPointSets.
		 */
		double [] hGeneral= {9,7.25,1.85,7.03,20.4,4.3};

		double [] h1= new double[hGeneral.length];


		for (int i=0;i<hGeneral.length;i++) {
			h1[i]=dim*hGeneral[i]/10;

		}


		/**
		 * For more details, refer to the class GenzfunctionTestWithBestPointSets.
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




		/**
		 * Number of simulations.
		 */
		 int nbreSimulateLMS=1000;

		 double []moyenne = new double[nbreSimulateLMS] ;
		 double [] variance = new double[nbreSimulateLMS];
		 double [] WafomResult=new double[nbreSimulateLMS];


		 /**
		  * Number of replicates for RQMC.
		  */
		 int replicates=20;
		 double []moyenneRQMC=new double[nbreSet];
		 double []varianceRQMC=new double[nbreSet];


		 for (int i=0;i<nbreSimulateLMS;i++) {
			 System.out.println("we are at simul number:"+(i+1));

			 /**
			  * We generate our point set.
			  */

			 DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
			 Sob.leftMatrixScramble(new MRG32k3a());

			 /**
			  * We calculate and store our WAFOM for each simulation.
			  */

			 Wafom waf=new Wafom (Sob,1,dim,w,k);

			 WafomResult[i]=waf.calcWafom();
			 PointSetRandomization rand = new RandomShift(new MRG32k3a());


			 for (int j=0;j<nbreSet;j++) {
				 stats1.init();






				 RQMCExperiment.simulReplicatesRQMC(new VarVsWafNaifSimul (dim,w, k,  nbreSet, aTotal[j], uTotal[j]),  Sob, rand, replicates,
						 stats1) ;

				 moyenneRQMC[j]=stats1.average();
				 varianceRQMC[j]=stats1.variance();


			 }

			 /**
			  * We calculate and store our variances for each simulation.
			  */
			 moyenne[i]=calculateMedian(moyenneRQMC);
			 variance[i]=calculateMedian(varianceRQMC);	

		 }



		 System.out.println("general results");
		 System.out.println("variance***********Wafom:");

		 for (int i=0;i<variance.length;i++) {
			 System.out.println(variance[i]+","+WafomResult[i]);

		 }



		 double [] rearrangeVariance=arrangeVecteur(variance);

		 double []rearrangeWafom = arrangeVecteur(WafomResult);
		 int [] indicesTri=getIndicesTri(WafomResult);


		 System.out.println("Here, we have sorted the best and worst Wafom and see the corresponding variances.");
		 System.out.println("Best: variance***********Wafom:");



		 int integerPart = (int) (nbreSimulateLMS/2);


		 for(int i=0;i<(integerPart);i++) {
			 System.out.println(variance[indicesTri[i]]+","+WafomResult[indicesTri[i]]);

		 }

		 System.out.println("plus mauvaises variance***********Wafom:");

		 for(int i=rearrangeWafom.length-1;(rearrangeWafom.length-(integerPart))<=i;i--) {
			 System.out.println(variance[indicesTri[i]]+","+WafomResult[indicesTri[i]]);

		 }



	}



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
				// Si la contrainte n'est pas respectée, reprenez la génération
				i--;
				continue;
			}

			aTotal[i][dim - 1] = aLast;
		}
		return aTotal;
	}





	public static double[] arrangeVecteur(double[] vecteur) {
		// Copie du vecteur original pour éviter de modifier l'original
		double[] vecteurCopie = Arrays.copyOf(vecteur, vecteur.length);

		// Tri du tableau accendant le 1er element est le plus petit le dernier et le plus grand
		Arrays.sort(vecteurCopie);

		return vecteurCopie;
	}
	public static int[] getIndicesTri(double[] vecteurOrigi) {
		int []indices=new int[vecteurOrigi.length];
		double [] vecteurOrdonne=arrangeVecteur(vecteurOrigi);
		//	double []vecteur1=Arrays.copyOf(vecteur, vecteur.length);
		for(int i=0;i<vecteurOrigi.length;i++) {
			for(int j=0;j<vecteurOrigi.length;j++) {
				if (vecteurOrdonne [i]==vecteurOrigi[j])
				{
					indices[i]=j;
				}
			}


		}
		return indices;

	}
}




