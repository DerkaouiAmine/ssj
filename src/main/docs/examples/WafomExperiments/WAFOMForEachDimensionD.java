package WafomExperiments;



import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;


/**
 * This method has not been outlined anywhere in my work. However, it's mentioned here just as an indicative note:
 * if we construct the best generator matrix for each dimension in terms of WAFOM and concatenate all the obtained matrices,
 * we do not obtain a significantly lower WAFOM than that obtained in the methods described in my work.
 */
public class WAFOMForEachDimensionD {
	 public static void main(String[] args) throws FileNotFoundException { 
	    	long startTime = System.currentTimeMillis();

		 int dim=3;
			int k=14;
	
		
			int w=30;
			
			int [][] matrixgenerat=new int[dim][k];
			
			
			
			double[]bestWAFOMS=new double[k*dim];
			double[]worstWAFOMS=new double[k*dim];
			
			int[] bestColonnes=new int [ k];
			int[] worstColonnes=new int [ k];

			int nbreTest=1000;

		        
		        
		     
				int[]Generamatrix = new int[k]; 
				
				
				/**
				 * The WAFOM is calculated for each dimension.
				 */

			    	for(int dimension=1;dimension<=dim;dimension++) {
			    		
			    		
			    	for(int d=0;d<k;d++) {
				double smallestWafom=(1<<k);
				double biggestWafom=-1;
			
			    	for(int i1=0;i1< nbreTest;i1++) {
			    		
			    		
			    	DigitalNetBase2 Sob = new SobolSequence(d+1,w,dimension);
				       Sob.leftMatrixScramble(new MRG32k3a());
				       /**
				        * We construct our optimal generator matrix for each dimension.
				        */

				       int [][] generatriOriginale=Sob.getGeneratorMatrices(d+1);
				       int element=generatriOriginale[dimension-1][d];
				       growingVector(Generamatrix, element, d);
				       int [] Generamatrix1 = Arrays.copyOfRange(Generamatrix, 0, d+1);
				       
				       DigitalNetBase2 Sob1 = new SobolSequence(d+1,w,1,Generamatrix1);//car je fais dimension par dimension donc 1 seul ce qui importe c'est les elemeents qui constitu la matrice generatrice

					      double [][]SobPoints = Sob1.formatPointsTab();
				     
					      Wafom lowWaf=new Wafom (Sob1,1,1,w,d+1);

					       double currentWafom =lowWaf.calcWafom1(SobPoints,d+1);
					       
					       /**
					        * We keep track of the best and worst WAFOM values.
					        */
					       
					       if (currentWafom < smallestWafom) {
					    	   smallestWafom = currentWafom;
	
					    	   Generamatrix[d]=element;
					    	   bestWAFOMS[d]=smallestWafom;
					    	   	bestColonnes[d]=element;
					    	   	//il faut d*index;

					       }

					       if (currentWafom > biggestWafom) {
					    	   biggestWafom= currentWafom;
	
	    
					    	   worstWAFOMS[d]=biggestWafom;
					    	   worstColonnes[d]=element;

					       }
				       

			    	}
			    	System.out.println(smallestWafom);
			    	
			    	
			    	}
			    	/**
			    	 * We store the best column for each dimension.
			    	 */

			    	matrixgenerat[dimension-1]=bestColonnes;
			    	System.out.println(Arrays.toString(bestColonnes));

			    	
	 }
			    	/**
			    	 * We construct our general matrix, which is the best of each dimension.
			    	 */

			    	int []vect=matrixToVector(matrixgenerat);
					
					DigitalNetBase2 Sob = new SobolSequence(k,w,dim,vect);
				  	 					
					 Wafom lowWaf=new Wafom (Sob,1,dim,w,k);
				        
				       double currentWafom =lowWaf.calcWafom();
				       
				       /**
				        * This WAFOM is significantly larger than all the methods studied in my work.
				        */

				       
				       System.out.println(currentWafom);
			  
	 
	 }
	 
	 
	 /**
	  * This method is essential for creating our vector for each dimension.
	  */

	 
	  public static void growingVector(int[] inputArray, int element, int place) {
	        int[] Vector = new int[inputArray.length];
	        inputArray[place]=element;
	        
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
	    
	 

	
}