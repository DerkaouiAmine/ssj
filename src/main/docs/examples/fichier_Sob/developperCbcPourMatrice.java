package fichier_Sob;



import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG31k3p;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;


//Methode2 ilfaut calculer le Waf pour chaque vecteur colonne et pas l'inserre a chaque fois pour avoir une matrice et comparerr les deux methodes

public class developperCbcPourMatrice {//marche tres mal ie meme si on ades wafom tres perit pour des 1D sa veut pas dire que la matrice aura un petit Wafom
	 public static void calc(int k) { 
	    	long startTime = System.currentTimeMillis();

		 int dim=6;
			//k=14 on peut allez jusqu a 7000 1.009
			//k=15 4.143 8heures pour 7000
			//k=16 19.340  pour 2000 c'est 11h
			//k=17 88.576 pour 1000 24h
			//k=18 420429 pour 1000 c'est 5jours
			//k= 19 1772094 m.secondes // pour 100 on a 2 jours
			//k=20 7413181 pour 100 c'est 8 jours
			
			
			
			//touver un moyen pour agir directement sur sob
			//ajouter une matrice a la declaration sobsequence (k,w,..,matricegenerat)

			
			int w=30;
			//je vais le modifier pour agir directement sur le tableau principale
			int nombdeDeColonnesExtraites =1 ;
			
			
			
			int nbreSimulateLMS=1000;
			int nbreSimulateDigi=100;
			
			double[]bestWAFOMS=new double[k*dim];
			double[]worstWAFOMS=new double[k*dim];
			int  bestColonne=0;
			int worstColonne=0;
			int[] bestColonnes=new int [ k];
			int[] worstColonnes=new int [ k];

			RandomStream randomStream = new MRG32k3a();

		        
		        
		        
				
				//int [][] matrix = new int[dim][k];
			int[][] Generamatrix= new int[dim][k];
			int[][] bestGeneramatrix= new int[dim][k];

			
			    //for (int l = 0; l < k; l++) {//l=colonne extraite
			    	for(int dimension=1;dimension<=dim;dimension++) {

			    	for(int d=0;d<k;d++) {
				double smallestWafom=(1<<k);
				double biggestWafom=-1;
				//int[] colonne1=new int[dim];
			    	for(int i1=0;i1<100;i1++) {
			    		//System.out.println("je suis a la "+i1+"de k="+(d+1)+"de dimesion="+dimension);
			    		
			    		
			    	DigitalNetBase2 Sob = new SobolSequence(d+1,w,dimension);
				       Sob.leftMatrixScramble(new MRG32k3a());
					
				       int [][] generatriOriginale=Sob.getGeneratorMatrices(d+1);
				  
				       
				       int element=generatriOriginale[dimension-1][d];
				       
				       
				       Generamatrix[dimension-1][d]=element;
				       
				       int [][] generat=extractSubMatrix(Generamatrix,dimension-1,d);
				       
				     /*  for(int []row:generat) {
				    	   System.out.println(Arrays.toString(row));
				       }*/
				       
				       int vect[]=matrixToVector( generat);
				       
				       
				       DigitalNetBase2 Sob1 = new SobolSequence(d+1,w,dimension,vect);

				     //sa plus vite ou la nouvelle que j'ai creer sur SobolSequqnce et comparer
				     			      double [][]SobPoints= Sob1.formatPointsTab();
				     			      
				     			     Wafom lowWaf=new Wafom (Sob1,1,dimension,w,d+1);
				 			        
								       double currentWafom =lowWaf.calcWafom1(SobPoints,d+1);
								       
								     if (currentWafom < smallestWafom) {
								           smallestWafom = currentWafom;
								           
								          // bestWAFOMS[d]=smallestWafom;
									    	bestColonne=element;
							 
								       }
								       
								       if (currentWafom > biggestWafom) {
								    	   biggestWafom= currentWafom;
								          // bestColonne= colonne1; // Mise à jour de la colonne correspondante
								           //System.out.println(i1+" fois de plus");
								           
								          // worstWAFOMS[d]=biggestWafom;
									    	worstColonne=element;
							 
								       }
				     			      
				       
				       
				       
				      
				       
				       
				     
			    	
	 }
			    	
			    	 Generamatrix[dimension-1][d]=bestColonne;
			  
	 
	 }
			    	}
			    	System.out.println("pour k="+k);
			    	for(int []vecteur:Generamatrix)
			    	{System.out.println(Arrays.toString(vecteur));}
			    	
			    	int []vect=matrixToVector(Generamatrix);
					
					DigitalNetBase2 Sob = new SobolSequence(k,w,dim,vect);
					
					
					 Wafom lowWaf=new Wafom (Sob,1,dim,w,k);
				        
				       double currentWafom =lowWaf.calcWafom();
				       
				       System.out.println(currentWafom);
			    	
	 }
	 
	 public static void main(String[] args) throws FileNotFoundException {
		 
		 for(int k=1;k<=23;k++)
		 {calc(k);}
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
	 
	 
	 
	 
	  public static void growingVector(int[] inputArray, int element, int place) {
	        int[] Vector = new int[inputArray.length];
	        inputArray[place]=element;
	        
	    }
	 
	 
	 

	   public static int extractElement(int[] vector, int index) {
	        if (index < 0 || index >= vector.length) {
	            throw new IndexOutOfBoundsException("Invalid index. Index must be between 0 and vector.length - 1.");
	        }

	        return vector[index];
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