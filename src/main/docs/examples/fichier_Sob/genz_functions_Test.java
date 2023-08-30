package fichier_Sob;
import java.util.Arrays;
import java.util.Random;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

public class genz_functions_Test {

    static final int NUM_SAMPLES = 10000000; // Augmentez le nombre d'échantillons
    static final int DIMENSION = 5; // Dimension de l'intégrale
    
    // Oscillatory function: cos(2*pi*u1 + sum(a[i] * x[i]))
    static double oscillatory(double[] x, double[] a,double [] u) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += a[i] * x[i];
        }
        return Math.cos(2 * Math.PI *u[0]  + sum);
    }
    
    // Product Peak function: prod(a[i]^-2 + (x[i] - u[i])^-2)
    static double productPeak(double[] x, double[] a, double[] u) {
        double product = 1.0;
        for (int i = 0; i < x.length; i++) {
        	double Result=1/(Math.pow(a[i], -2)+Math.pow(x[i]-u[i],2));
        	
            product *= Result;
        }
        return product;
    }

    // Corner Peak function: (1 + sum(a[i] * x[i]))^(-s + 1)
    static double cornerPeak(double[] x, double[] a) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += a[i] * x[i];
        }
        return Math.pow(1 + sum, -x.length + 1);
    }

    // Continuous function: exp(-sum(a[i] * |x[i] - u[i]|))
    static double continuous(double[] x, double[] a,double [] u) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += a[i] * Math.abs(x[i] - u[i]);
        }
        return Math.exp(-sum);
    }

    // Discontinuous function: 0 if x[i] > u[i], otherwise exp(sum(a[i] * x[i]))
    static double discontinuous(double[] x, double [] a,double[] u) {
        for (int i = 0; i < x.length; i++) {
            if ((x[1] > u[1]) ||(x[2] > u[2])) {
                return 0.0;
            }
        }
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += a[i] * x[i];
        }
        return Math.exp(sum);
    }
    

    static double gaussian(double[] x, double[] a, double[] u) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += a[i] * a[i] * Math.pow(x[i] - u[i], 2);
        }
        return Math.exp(-sum);
    }

    public static void main(String[] args) {
    	
	    RandomStream random = new MRG32k3a();
	    int nbreSet = 20;
	    
	   
       
        double[][] uTotal = new double[nbreSet][DIMENSION]; // Constants u_i
        
        double [] integralValue=new double[nbreSet]; 
        
        
        //hj pour s=10 ,OSc=h0,productPEak=h1,CornerPeak=h2,Gaussi=h3,Continious=h4,Discontin=h5
        double [] hGeneral= {9,7.25,1.85,7.03,20.4,4.3};//  voir reference  pour s=10(Barthelmann.NOVAK)

        
        double [] h1= new double[hGeneral.length];//3/10*ce que a donne pour s=10
        		
        		
        for (int i=0;i<hGeneral.length;i++) {
        	h1[i]=DIMENSION*hGeneral[i]/10;

        }
       
        
        //choisir le rang de la fonction choisie
        double h = h1[3];
       
       for (int i=0;i<nbreSet;i++) {
        random.nextArrayOfDouble(uTotal[i], 0, DIMENSION);
        }
    
       

       
      double [][] aTotal= generateSetForA (nbreSet, h, DIMENSION);// Constants a_i
 
 
      
      for(double [] tab : aTotal) {
    	  System.out.println(Arrays.toString(tab));
    	  
      }
      System.out.println("les U");
      for(double [] tab1 : uTotal) {
    	  System.out.println(Arrays.toString(tab1));
    	  
      }
      
  
        
 System.out.println("POUR I(f)");
        
      for (int i1=0;i1<nbreSet;i1++) {
//System.out.println("je suis a la"+(i1+1));
        double integralSum = 0.0;

        for (int i = 0; i < NUM_SAMPLES; i++) {
            double[] x = new double[DIMENSION];
            for (int j = 0; j < DIMENSION; j++) {
                x[j] = random.nextDouble();
            }
           integralSum += gaussian(x, aTotal[i1], uTotal[i1]);
            // integralSum +=  oscillatory(x, aTotal[i1], uTotal[i1]);

        }

        double volume = Math.pow(1.0, DIMENSION); // Volume de l'hypercube unité
        integralValue[i1] = (integralSum / NUM_SAMPLES) * volume;
        
      }
      
      
     System.out.println("VAleurs de l'integrale");
    	  System.out.println(Arrays.toString(integralValue));
      
    	  double a=calculateMedian(integralValue);
    	  System.out.println("LA mediane est :"+a);
    	  
    	  
    	  
    	  
    	  
    	  
    	  
    	  
    	  
    	  
    	  
    	  // MAINTENANT POUR LE CAS WAFOM
    	  
    	  
    	  System.out.println("POUR IN(f)");
          
    	  
    	  
    		int dim=3;
    	    int dimension=0;
    		int k=14;//k=23 est le maximum dans notre cas
    		int w=30;
    		
    		double cste=Math.pow(2, -w-1);
    		
    		double [] mediane=new double[k];
    		
    		Tally stats=new Tally("Simul");
    		
    		double [] [] ResultatsIn=new double[k][nbreSet];
    	  
    	  //MEILLEUR 7000
    		/* int[][] Generamatrix =  {
    				 {1046994105,398307896,181181815,95054231,34220517,21484503,12394050,4574640,2364507,1049891,926899,518558,221597,94087,34218,30883,14575,4451,2382,1950,554},
    				 {1069297368,735226078,767106076,843169012,899961314,941641228,602368976,667437869,788484991,955352261,742762285,611944903,919851144,779148930,809466920,789232013,584593252,1051409679,983178797,677711189,996379055},
    				 {1065712118,528293238,620731993,863648038,471982456,797987891,577453692,437494779,1031238497,885621485,347984506,622078176,850186234,393980218,972996136,1044507081,465061887,755970565,776187327,289374307,845703665}
    		 };*/
    		
    		//meilleur 100000 pour s=5
    		int[][] Generamatrix =  {
    		 {1050571877, 392693088, 210662097, 117269155, 49184842, 32843177, 15551431, 6506393, 2992302, 1435308, 715399, 391043, 255604, 105156},
             {1060625895, 585690982, 615590625, 1044883373, 798396084, 953062521, 927185581, 925485110, 1006437641, 1029032659, 946002992, 759664060, 864799114, 668596113},
             {1028833161, 513977673, 611061055, 859066087, 511182609, 841135827, 796064304, 410216384, 629229509, 713180353, 472177215, 651517590, 786105168, 409650319 },
             {1056559733, 867065376, 598896813, 454886611, 213422505, 871874645, 451894096, 733159149, 589424565, 808946694, 428404476, 143783704, 816636780, 492612364},
             {1011151336, 401782358, 586968320, 217427154, 455800134, 991286047, 987396119, 580240759, 335152366, 660700410, 224697316, 385137884, 839855277, 564125023}
         };
    		 
    		 
    		 
    		 
    		  double [] integralValue1=new double[nbreSet]; 

    		 
    		 for (int k1=1;k1<=Generamatrix[0].length;k1++) {
    		      
    			 
    			 //System.out.println("cas "+(k1));
    		        
    			 stats.init();
    		
    			 int [][] generat=extractSubMatrix(Generamatrix,Generamatrix.length-1,k1-1);
    			 
    			 
    			 int []vect=matrixToVector(generat);
    			 
    			 DigitalNetBase2 Sob = new SobolSequence(k1,w,dim,vect);
    			// for (int m=0;m<100;m++) {
					Sob.addRandomShift(new MRG32k3a());
    		
    			 double [][] pointSet=Sob.formatPointsTab();
    			 
    			 

					int nmbrePoints=(1<<k1);
    			 
    		      for (int i1=0;i1<nbreSet;i1++) {
    		    	          double integralSum = 0.0;
    		    	          
		    	        		

		    					
		    					

    		    	          for (int i = 0; i < nmbrePoints; i++) {
    		    	        	  
    		    				
    		    	        	  
    		    	        	  
    		    	        	  
    		    	        	  
    		    	              double[] x = pointSet[i];
    		    	              
    		    	            /*  for (int l = 0; l < x.length; l++) {
    		    	                  x[l] += cste;
    		    	              }*/

    		    	              
    		    	          
    		    	              integralSum += gaussian(x, aTotal[i1], uTotal[i1]);
    		    	              
    		    	             // integralSum +=   oscillatory(x, aTotal[i1], uTotal[i1]);
    		    	          }
    		    	          double volume = Math.pow(1.0, DIMENSION); // Volume de l'hypercube unité

    		    	          integralValue1[i1] = (integralSum / nmbrePoints) * volume;
    		    	          ResultatsIn[k1-1][i1] = (integralSum / nmbrePoints) * volume;

    		    	          
    		    	        }
    		       System.out.println("VAleurs de l'integrale pour k= "+(k1));
    	      	  System.out.println(Arrays.toString(integralValue1));
    	        
    	      	  mediane[k1-1]=calculateMedian(integralValue1);
    	      	  System.out.println("LA mediane est :"+mediane[k1-1]);
    	      	  
    		 }
    		    	        
    		    	        
    		    	      
    		    	      	  
    			 
    			 
    			 
    			 
    			 
    			 
    			 
    		
    		
    
    		// la mediane est calculer a chaque etape
    		 
    		 System.out.println("LES VALEURS A TRACER");
    		 for (int i2=0;i2<mediane.length;i2++) {

    			 double a1=Math.abs(a-mediane[i2])/Math.abs(a);
    			 //double a1=Math.abs(mediane[i2]);
        		 System.out.println(Math.log10(a1));

    			 
    		 }
    		 
    		 //la mediane apres avoir calculer le tt
  		 /*double [] aTracer= new double[k];

    		 
    		 for (int m=0;m<k;m++) {
        		 double [] Resultats= new double[k];

    		 
    		 for (int l=0;l<nbreSet;l++) {
    			 
    			 Resultats[l]=(integralValue[l]-ResultatsIn[m][l])/integralValue[l];
    			 
    		 }
    		 aTracer[m]=calculateMedian(Resultats);
    		 }
    		 
    		 for (int i=0;i<aTracer.length;i++) {
    			 System.out.println(Math.log10(aTracer[i]));
    		 }*/
    		 
    		 
    		 
  
    		 
    		 
    		 
        
        
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
