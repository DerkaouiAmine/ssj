package fichier_Sob;

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
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;






//reste a teste le biais pour a variance  on prend encore lms +100 digital shift et encore 100 digital shift on vois si sa change ou pas pour voir si y'
public class randomShiftAvecUnSeulWafom {
	 public static void main(String[] args) throws FileNotFoundException  { 
		 long startTime = System.currentTimeMillis();

	int dim=3;
    int dimension=0;
	int k=16;//k=23 est le maximum dans notre cas
	int w=30;
	
	
	 double functionValue=0;
	 

	
	Tally stats=new Tally("Simul");
	Tally stats1=new Tally("Vari");
	
	int nbreSimulateDigi=1000;

	double []moyenne = new double[nbreSimulateDigi] ;
	double [] variance = new double[nbreSimulateDigi];
	

	double [][] pointSet=new double[1<<k][dim];
	
	
	
	
		
		pointSet=new double[1<<k][dim];
		//la meilleur du CBC 7000
		 int[][] Generamatrix =  {
				 {1046994105,398307896,181181815,95054231,34220517,21484503,12394050,4574640,2364507,1049891,926899,518558,221597,94087,34218,30883,14575,4451,2382,1950,554},
				 {1069297368,735226078,767106076,843169012,899961314,941641228,602368976,667437869,788484991,955352261,742762285,611944903,919851144,779148930,809466920,789232013,584593252,1051409679,983178797,677711189,996379055},
				 {1065712118,528293238,620731993,863648038,471982456,797987891,577453692,437494779,1031238497,885621485,347984506,622078176,850186234,393980218,972996136,1044507081,465061887,755970565,776187327,289374307,845703665}
		 };
		
		/*// Naive sob+LMS 1 seul fois
		int[][] Generamatrix =  {
				{848380610, 335182282, 224252740, 114343632, 65672281, 24742595, 12226154, 4881155, 2344009, 1392392, 550815, 438448, 217232, 122380, 53520, 17720},
				{588942007, 1062796041, 680664071, 867509452, 560651987, 1015661661, 711420154, 816648793, 589766700, 1065211686, 681736757, 865850489, 559281240, 1018122941, 710232155, 813750047},
				{854021773, 397020222, 745672510, 1068526262, 502085461, 699714164, 871020778, 363658652, 792682128, 1040266498, 532157002, 707019989, 851723305, 395568882, 744898186, 1066738709}
	 };*/
		 
	/*	 //best matrices element by elements
		 int[][] Generamatrix =  {
{1073607198, 526206411, 268387157, 131651045, 66023764, 31447671, 16481707, 6688102, 4119781, 1531225, 714474, 441267, 217757, 104204, 44252, 19134},
{1072418975, 550987606, 1041731949, 768220889, 887700036, 1000103701, 777000097, 668099810, 934009139, 1020571374, 959400069, 591771954, 890317156, 872020837, 1004686711, 983487101},
{1066928913, 535688917, 550159880, 975246905, 369097981, 1017094843, 834401849, 530987773, 646108930, 608526229, 515278874, 941802869, 839001281, 372432652, 897826448, 666870939}

		 };*/
		int []vect=matrixToVector(Generamatrix);
		
		
		
	
		
		
		
			for (int j=0;j<nbreSimulateDigi;j++) {
				stats.init();
				System.out.println("Je suis a la "+j);

				DigitalNetBase2 Sob = new SobolSequence(k,w,dim,vect);
				functionValue=0;
		    	  PointSetRandomization rand1 = new RandomShift(new MRG32k3a());
		    	  rand1.randomize(Sob);
				//Sob.addRandomShift(new MRG32k3a());
				
				
				pointSet=Sob.formatPointsTab();
			
				for (int i1 = 0; i1 < pointSet.length; i1++) {
				    double somme = 0; // Déplacer l'initialisation ici, à l'intérieur de la boucle externe
				    RandomStream randomStream1 = new MRG32k3a();
				    for (int i2 = 0; i2 < pointSet[0].length; i2++) {
				        somme += Math.pow(randomStream1.nextDouble(), 2)* Math.pow((pointSet[i1][i2] - randomStream1.nextDouble()),2);
				    }
				    functionValue = Math.exp(-somme);
				    stats.add(functionValue);
				}
				moyenne[j]=stats.average();
			    variance[j]=stats.variance();
				
	}
				
		 //   System.out.println("je suis a la "+i);
			
	
	
	System.out.println("*************MOYENNE********************");
	for(int i=0;i<moyenne.length;i++)
	{System.out.println(moyenne[i]);}
	System.out.println("*************Variance********************");
	for(int i=0;i<variance.length;i++)
	{System.out.println(variance[i]);
	stats1.add(variance[i]);}
	System.out.println("*********************************");
	
	
	
	int [] indicesTri=getIndicesTri(variance);

	System.out.println("Best variance ");

	//System.out.println("Les 20 meilleur Wafom:");
	//je ppuvais utiliser juste les indices
	for(int i=0;i<10;i++) {
		//System.out.println(rearrangeWafom[i] +" variance correspanante : "+variance[indicesTri[i]]);
		System.out.println(variance[indicesTri[i]]);

	}
	
	System.out.println("Worst variance ");
	for(int i=indicesTri.length-1;i>(indicesTri.length-10);i--) {
		//System.out.println(rearrangeWafom[i] +" variance correspanante : "+variance[indicesTri[i]]);
		System.out.println(variance[indicesTri[i]]);

	}
	
	System.out.println(stats.report());

	
	 


 }

	 
	 public static double[] ExtractVecteur(double[][] resultatMatr,int j)
		
		{double[] vj = new double[resultatMatr.length];
			for (int i = 0; i < resultatMatr.length; i++) {
	            vj[i] = resultatMatr[i][j];
	        }
			return vj;
		}
		
	 
	 
	 
		public static double evaluateFunction(double[] point) {
	        
	    	RandomStream randomStream = new MRG32k3a();
	    	double Somme=0;
	    	for (int i=0;i<point.length;i++)
	    	{
	    		Somme+=(1 + randomStream.nextDouble() * (point[i] - 0.5));
	    	}
	        return Somme/point.length;
		
		}
		
	    public static double calculateNorm(double[] vector) {
	        double sumOfSquares = 0.0;
	        
	        for (double component : vector) {
	            sumOfSquares += Math.pow(component, 2);
	        }
	        
	        return Math.sqrt(sumOfSquares);
	    }
	    
	    public static int[] findMinMax(double [] v)
	    {
	        int minIndex = 0;
	        int maxIndex = 0;
	        double minValue = v[0];
	        double maxValue = v[0];
	        
	        for (int i = 1; i < v.length; i++) {
	            if (v[i] < minValue) {
	                minValue = v[i];
	                minIndex = i;
	            }
	            
	            if (v[i] > maxValue) {
	                maxValue = v[i];
	                maxIndex = i;
	            }
	        }
	        int[] indices = {minIndex, maxIndex};
	        return indices;
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
	
	
	 

