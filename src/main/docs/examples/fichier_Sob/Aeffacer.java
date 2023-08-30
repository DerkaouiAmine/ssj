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
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;






//reste a teste le biais pour a variance  on prend encore lms +100 digital shift et encore 100 digital shift on vois si sa change ou pas pour voir si y'
public class Aeffacer implements MonteCarloModelDouble {
	static int dim;
	double sum;
	 int nbreSet;
	 int w;
	 static double []a;
	 int k;
	 static double []u;

	// Constructor.
	public Aeffacer (int dim,int w,int k, int nbreSet,double[] a, double[] u) {
		this.dim = dim; 	this.nbreSet= nbreSet;
		this.w=w;			this.k=k;
		this.a=a;			this.u=u;
	}

	// Generates and returns X, without IS.
	public void simulate (RandomStream stream) {
		sum= 0.0;
		double [] x=new double [dim];
		for (int j = 0; j < dim; j++) {
			x[j]= stream.nextDouble();
			
		}
		
		sum=genz_functions.gaussian(x, a, u);
		//sum=genz_functions.oscillatory(x, a, u);//changer le h
	}
	
	
	

	
	

	// Generates and returns X, without IS.
	public double getPerformance () {
		return sum;
	}

	// Descriptor of this model.
	public String toString () {
		return "Test function for MC and RQMC: product of exponentials and cosine functions.";
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
		        
		        double [] integralValue=new double[nbreSet]; 
		        
		        
		        //hj pour s=10 ,OSc=h0,productPEak=h1,CornerPeak=h2,Gaussi=h3,Continious=h4,Discontin=h5
		        double [] hGeneral= {9,7.25,1.85,7.03,20.4,4.3};//  voir reference  pour s=10(Barthelmann.NOVAK)

		        
		        double [] h1= new double[hGeneral.length];//3/10*ce que a donne pour s=10
		        		
		        		
		        for (int i=0;i<hGeneral.length;i++) {
		        	h1[i]=dim*hGeneral[i]/10;

		        }
		       
		        
		        //choisir le rang de la fonction choisie
		        double hi = h1[3];
		       
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
		      int nbreSimulateLMS=100;
		      int nbreSimulateDigi=100;
		      double []moyenne = new double[nbreSimulateLMS] ;
		  	double [] variance = new double[nbreSimulateLMS];
		  	double [] WafomResult=new double[nbreSimulateLMS];
		  	int replicates=20;
		  	double []moyenneRQMC=new double[nbreSet];
			double []varianceRQMC=new double[nbreSet];
		  	
	

		
		 
		WafomsStorage waf=new	WafomsStorage(dim, 1,1);
		int [][]generat=waf.getWafomsTabl();
		
		
		int [] vect=matrixToVector(generat);
		DigitalNetBase2 Sob = new SobolSequence(k,w,dim,vect);
		System.out.println("sa ne veux rien dire car on a toujours le meme wafom on peut tester grace au k ");

   	  PointSetRandomization rand = new RandomShift(new MRG32k3a());
   	  for (int m=0;m<nbreSimulateDigi;m++) {
System.out.println("sans mediane pour j="+m);
		
		  for (int j=0;j<nbreSet;j++) {
			  stats1.init();

	    	  
	    	  
	    	  
			
	    	  
	    	  RQMCExperiment.simulReplicatesRQMC(new MCRQMCTEST (dim,w, k,  nbreSet, aTotal[j], uTotal[j]),  Sob, rand, replicates,
						stats1) ;
	      
	      moyenneRQMC[j]=stats1.average();
	      varianceRQMC[j]=stats1.variance();
	      //System.out.println("Moyenne = "+stats1.average()+" VAriance = "+stats1.variance());
	      
System.out.println(moyenneRQMC[j]+","+varianceRQMC[j]);
		
	      }
		  moyenne[m]=calculateMedian(moyenneRQMC);
		    variance[m]=calculateMedian(varianceRQMC);	
	
   	  }
	
	
System.out.println("mediane des variance toujours pour le memes Wafom");
	 
	System.out.println(Arrays.toString(variance));
	
	
	System.out.println("mediane des moyenne toujours pour le memes Wafom");
	 
	System.out.println(Arrays.toString(moyenne));
	
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
	    
	    
	    
	/*public static double [][] generateSetForA (int numSets, double targetSum, int dim) {
		double[][] aTotal = new double[numSets][dim];
	RandomStream random = new MRG32k3a();

	for (int i = 0; i < numSets; i++) {
	    double a1 = random.nextDouble();
	    double a2 = random.nextDouble();
	    double a3=random.nextDouble();
	    double a4=random.nextDouble();

	    double a5 = targetSum - a1 - a2 -a3-a4;

	    if (a5 < 0 || a5 > 1) {
	        // Si la contrainte n'est pas respectée, reprenez la génération
	        i--;
	        continue;
	    }
	    // Ajout de [a1, a2, a3] à aa[i]
	    double[] newSet = { a1, a2, a3 ,a4,a5};
	    System.arraycopy(newSet, 0, aTotal[i], 0, newSet.length);
	    
	   // System.out.println("Set " + (i + 1) + ": a1 = " + a1 + ", a2 = " + a2 + ", a3 = " + a3);
	}
	return aTotal;
	}*/
	    
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
	    
	    }
	
	
	 

