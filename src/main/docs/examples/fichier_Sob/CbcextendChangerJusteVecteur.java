package fichier_Sob;


// ici la difference est que a chaque fois on garde pas l'ancienne on test juste avec la colonne  
import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;


//Methode2 ilfaut calculer le Waf pour chaque vecteur colonne et pas l'inserre a chaque fois pour avoir une matrice et comparerr les deux methodes

public class CbcextendChangerJusteVecteur {//je prend par colonne
	 public static void main(String[] args) throws FileNotFoundException { 
	    	long startTime = System.currentTimeMillis();

		 int dim=3;
		    int dimension=0;
			int k=23;
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
			
			double[]bestWAFOMS=new double[k];
			double[]worstWAFOMS=new double[k];
			int [] bestColonne=new int[dim];
			
			int[][] bestColonnes=new int [ k][dim];
			int[][] worstColonnes=new int [ k][dim];


		        
		        
		        
				
				//int [][] matrix = new int[dim][k];
				int[][] Generamatrix = new int[dim][0]; // Matrice de sortie avec zéro colonne initialement

			    //for (int l = 0; l < k; l++) {//l=colonne extraite
			    	
			    	for(int d=0;d<k;d++) {
			    		
			    		
			    		System.out.println("Je calcule pour k="+(d+1));

				double smallestWafom=(1<<k);
				double biggestWafom=-1;
				//int[] colonne1=new int[dim];
			    	for(int i1=0;i1<1000;i1++) {
			    	DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
				       Sob.leftMatrixScramble(new MRG32k3a());
					
					int [][]	tableauPrincipalGenerLMS =  Sob.getGeneratorMatrices(k);//ici sa sera un i au lieu de 0

			    	
					 Generamatrix=tableauPrincipalGenerLMS;
					//int[] colonne1=new int[dim];

					int[] colonne1 = extraireColonne(tableauPrincipalGenerLMS, d);//ici i1 a la place de 0
			     
			    if (i1==0) {
			        Generamatrix = concatenateMatrixAndColumn(Generamatrix, colonne1);}
			        else {
				           Generamatrix=replaceColumn(Generamatrix,colonne1, d);
			        }

			    
			  int vect[]=matrixToVector(Generamatrix);

				DigitalNetBase2 Sob1 = new SobolSequence(d+1,w,dim,vect);

//sa plus vite ou la nouvelle que j'ai creer sur SobolSequqnce et comparer
			      double [][]SobPoints= Sob1.formatPointsTab();
			      
			      
		
			        
			      
			     
			       
			      //  double[][] SobPoints1=Sob.formatPointsTab();

			     //   System.out.println("SOBOL: "+ SobPoints.length+"GENERA LIGNE: "+Generamatrix.length+"GENERA Colo: "+Generamatrix[0].length);

			/*  System.out.println("Sobol's Evolution1**********");
			        for (int i = 0; i <  SobPoints.length; i++) {
			            for (int j = 0; j <  SobPoints[0].length; j++) {
			                System.out.print( SobPoints[i][j] + " ");
			            }
			            System.out.println();
			        }
			        */
				   
				     
			      
			        Wafom lowWaf=new Wafom (Sob1,1,dim,w,d+1);
			        
				       double currentWafom =lowWaf.calcWafom1(SobPoints,d+1);
				       
				       if (currentWafom < smallestWafom) {
				           smallestWafom = currentWafom;
				          // bestColonne= colonne1; // Mise à jour de la colonne correspondante
				           //System.out.println(i1+" fois de plus");
				           
				           bestWAFOMS[d]=smallestWafom;
					    	bestColonnes[d]=colonne1;
			 
				       }
				       
				       if (currentWafom > biggestWafom) {
				    	   biggestWafom= currentWafom;
				          // bestColonne= colonne1; // Mise à jour de la colonne correspondante
				           //System.out.println(i1+" fois de plus");
				           
				           worstWAFOMS[d]=biggestWafom;
					    	worstColonnes[d]=colonne1;
			 
				       }
//System.out.println("FINI");
			        
				    /*	System.out.println("best WAFOMS");
				    	System.out.println(Arrays.toString( bestWAFOMS));
				    	
				    	
				        System.out.println("BEST COLOMNS");
				        for (int i = 0; i <  bestColonnes[0].length; i++) {
				            for (int j = 0; j <  bestColonnes.length; j++) {
				                System.out.print(bestColonnes[j][i] + " ");
				            }
				            System.out.println();
				        }
			    	
				        
				    	System.out.println("worst WAFOMS");
				    	System.out.println(Arrays.toString( worstWAFOMS));
				    	
				    	
				        System.out.println("worst COLOMNS");
				        for (int i = 0; i <  worstColonnes[0].length; i++) {
				            for (int j = 0; j <  worstColonnes.length; j++) {
				                System.out.print(worstColonnes[j][i] + " ");
				            }
				            System.out.println();
				        }
			    	
			    	}
			    	*/
			    	
			    	
			    	
			    	
			    	
			    
			    	
			 
	 
		  
			        
			    	//long endTime = System.currentTimeMillis();
			    //	long executionTime = endTime - startTime;
			    	
			    	//System.out.println(executionTime);
			    	

}
			    	
			    	
			    	
			   System.out.println("best WAFOMS");
			    	System.out.println(Arrays.toString( bestWAFOMS));
			    	
			    	
			        System.out.println("BEST COLOMNS");
			        for (int i = 0; i <  bestColonnes[0].length; i++) {
			            for (int j = 0; j <  bestColonnes.length; j++) {
			                System.out.print(bestColonnes[j][i] + " ");
			            }
			            System.out.println();
			        }
		    	
			        
			    	System.out.println("worst WAFOMS");
			    	System.out.println(Arrays.toString( worstWAFOMS));
			    	
			    	
			        System.out.println("worst COLOMNS");
			        for (int i = 0; i <  worstColonnes[0].length; i++) {
			            for (int j = 0; j <  worstColonnes.length; j++) {
			                System.out.print(worstColonnes[j][i] + " ");
			            }
			            System.out.println();
			        }
		    	
		    	}
			    	
			    	
			    	
			    	}
			    	
	 
	 public static int[][] concatenateMatrixAndColumn(int[][] matrix, int[] colonne1) {
		    int numRows = matrix.length;
		    int numCols = matrix[0].length;

		    if (numRows != colonne1.length) {
		        throw new IllegalArgumentException("Le nombre de lignes de la matrice doit être égal à la taille du vecteur colonne.");
		    }

		    int[][] resultMatrix = new int[numRows][numCols + 1];

		    for (int i = 0; i < numRows; i++) {
		        for (int j = 0; j < numCols; j++) {
		            resultMatrix[i][j] = matrix[i][j];
		        }
		        resultMatrix[i][numCols] = colonne1[i];
		    }

		    return resultMatrix;
		}
	 
	 public static int[] extraireColonne(int[][] tableau, int j) {
		    int nbLignes = tableau.length;
		    int[] colonne = new int[nbLignes];

		    for (int i = 0; i < nbLignes; i++) {
		        colonne[i] = tableau[i][j];
		    }

		    return colonne;
		}

	 
	 
	    public static int[][] extraireColonnes(int[][] tableau, int colonnesExtraites) {
	        int[][] colonnes = new int[tableau.length][colonnesExtraites];

	        for (int i = 0; i < tableau.length; i++) {
	            for (int j = 0; j < colonnesExtraites; j++) {
	                colonnes[i][j] = tableau[i][j];
	            }
	        }

	        return colonnes;
	    }
	    
	    
	    public static int[][] replaceColumn(int[][] matrix, int[] newColumn, int columnIndex) {
	        if (columnIndex < 0 || columnIndex >= matrix[0].length) {
	            throw new IllegalArgumentException("L'index de colonne spécifié est invalide.");
	        }

	        if (newColumn.length != matrix.length) {
	            throw new IllegalArgumentException("Le nouveau vecteur colonne doit avoir la même taille que le nombre de lignes de la matrice.");
	        }

	        for (int i = 0; i < matrix.length; i++) {
	            matrix[i][columnIndex] = newColumn[i];
	        }

	        return matrix;
	    }
	    
	    
	    public static int[][] vectorToMatrix(int[] vector) {
	        int numRows = vector.length;
	        int[][] matrix = new int[numRows][1];

	        for (int i = 0; i < numRows; i++) {
	            matrix[i][0] = vector[i];
	        }

	        return matrix;
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