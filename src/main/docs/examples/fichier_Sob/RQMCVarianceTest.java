package fichier_Sob;

import java.io.*;
import java.util.Arrays;

import umontreal.ssj.hups.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.mcqmctools.*;



public class RQMCVarianceTest implements MonteCarloModelDouble {
	static int dim;
	double sum;
	 int nbreSet;
	 int w;
	 static double []a;
	 int k;
	 static double []u;

	// Constructor.
	public RQMCVarianceTest (int dim,int w,int k, int nbreSet,double[] a, double[] u) {
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

	public static void main(String[] args) throws IOException {
		int k=21;
		int w=30;
		int dim=3;
		
		Tally stats=new Tally("Simul MC");
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
	 
	      
	
	    //pour RQMC+WAFOM
	      
	      
	 		/* RandomStream random2 = new MRG32k3a();
	 		PointSetRandomization rand = new RandomShift(random2);*/
	      
	      System.out.println("******************Methode RQMC****************");
	      
	      int s = dim;
	       // int methode=6;
	        
	        for (int methode=0;methode<7;methode++) {
	  	      System.out.println("******************Methode numero"+methode+"****************");
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
	      /* for (int[] l:Generamatrix ) {
	        	System.out.println(Arrays.toString(l));
	        }*/

	    /*  int[][] Generamatrix =  {
				 {1046994105,398307896,181181815,95054231,34220517,21484503,12394050,4574640,2364507,1049891,926899,518558,221597,94087,34218,30883,14575,4451,2382,1950,554},
				 {1069297368,735226078,767106076,843169012,899961314,941641228,602368976,667437869,788484991,955352261,742762285,611944903,919851144,779148930,809466920,789232013,584593252,1051409679,983178797,677711189,996379055},
				 {1065712118,528293238,620731993,863648038,471982456,797987891,577453692,437494779,1031238497,885621485,347984506,622078176,850186234,393980218,972996136,1044507081,465061887,755970565,776187327,289374307,845703665}
		 };*/
	 		 
			Tally stats1=new Tally("Simul RQMC");
			
			
			
			//a quoi sa sert replicates???
			int replicates=100;
			
			double []moyenneRQMC=new double[nbreSet];
			double []varianceRQMC=new double[nbreSet];
			double [] medianeDesVariance=new double[Generamatrix[0].length];
			
			double [][] moyennesRQMCGenerale=new double[Generamatrix[0].length][nbreSet];
			double [][] varianceRQMCGenerale=new double[Generamatrix[0].length][nbreSet];

			
			 for (int k1=1;k1<=Generamatrix[0].length;k1++) {
			     // System.out.println("--------------------------pour k= "+k1);

				 
				 
				 
				 int [][] generat=extractSubMatrix(Generamatrix,Generamatrix.length-1,k1-1);
  			 
  			 
  			 int []vect=matrixToVector(generat);
  			 
  			 DigitalNetBase2 Sob = new SobolSequence(k1,w,dim,vect);
					//Sob.addRandomShift(new MRG32k3a());

			
		     for (int i=0;i<nbreSet;i++) {
			    //  System.out.println("-----pour le set numero ="+(i+1));

		    	  stats1.init();
		    	  
		    	  PointSetRandomization rand = new RandomShift(new MRG32k3a());
		    	  
					// DigitalNetBase2 p = new SobolSequence(k,w,dim);

		    	  
		    	  RQMCExperiment.simulReplicatesRQMC(new MCRQMCTEST (dim,w, k,  nbreSet, aTotal[i], uTotal[i]),  Sob, rand, replicates,
							stats1) ;
		      
		      moyenneRQMC[i]=stats1.average();
		      varianceRQMC[i]=stats1.variance();
		     // System.out.println("Moyenne = "+stats1.average()+" VAriance = "+stats1.variance());

			
		      }
		     moyennesRQMCGenerale[k1-1]=moyenneRQMC;
		     varianceRQMCGenerale[k1-1]=varianceRQMC;
		     
		     medianeDesVariance[k1-1]=calculateMedian(varianceRQMC);
		     System.out.println( medianeDesVariance[k1-1]);

			
			 }
			 BestOrWorstWaf++;
	        }
	        
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



}
