package fichier_Sob;

import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;


//Methode2 ilfaut calculer le Waf pour chaque vecteur colonne et pas l'inserre a chaque fois pour avoir une matrice et comparerr les deux methodes

public class CBCtresLent{
	 public static void main(String[] args) throws FileNotFoundException { 
	    	long startTime = System.currentTimeMillis();

		 int dim=3;
		    int dimension=0;
			int k=18;
			int w=30;
			//je vais le modifier pour agir directement sur le tableau principale
			int nombdeDeColonnesExtraites =1 ;
			
			
			
			int nbreSimulateLMS=1000;
			int nbreSimulateDigi=100;
			
			double[]WAFOMS=new double[k];
			
			int [] bestColonne=new int[dim];
			
			int[][] bestColonnes=new int [ k][dim];
			

		        
		        
		        
				
				//int [][] matrix = new int[dim][k];
				int[][] Generamatrix = new int[dim][0]; // Matrice de sortie avec zéro colonne initialement

			    //for (int l = 0; l < k; l++) {//l=colonne extraite
			    	
			    	for(int d=0;d<k;d++) {
				double smallestWafom=(1<<3);
				int[] colonne1=new int[dim];
			    	for(int i1=0;i1<1;i1++) {
			    	DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
				       Sob.leftMatrixScramble(new MRG32k3a());
					
					//int[][] []tableauPrincipalGenerLMS = new int[1][][];

					//tableauPrincipalGenerLMS[0] =  Sob.getGeneratorMatrices(k);//ici sa sera un i au lieu de 0
					int [][]	tableauPrincipalGenerLMS =  Sob.getGeneratorMatrices(k);//ici sa sera un i au lieu de 0

			    	
					///faire ATTENTION LE K DANS 			        double [][]SobPoints= Sob.getSobolPointSetForGeneratMatrix(1, w, Generamatrix );//ici k=1
                       //DOIT ETRE LE MEME QUE LE NOMBRE DE COLONNE DE LA GENERATR MATRIX
					
					
					
			    	
			     colonne1 = extraireColonne(tableauPrincipalGenerLMS, d);//ici i1 a la place de 0
			    if (i1==0) {
			        Generamatrix = concatenateMatrixAndColumn(Generamatrix, colonne1);}
			        else {
				           Generamatrix=replaceColumn(Generamatrix,colonne1, d);

			       }
			       
			        
			        
			        
			        
			       /* if (i1==10000-1) {
					System.out.println("BBBBBBBBBBBBBBBBBBBBBBBBBBBB");

			        for (int i = 0; i <  Generamatrix.length; i++) {
			            for (int j = 0; j <  Generamatrix[0].length; j++) {
			                System.out.print( Generamatrix[i][j] + " ");
			            }
			            System.out.println();
			        }
			        
					System.out.println("BBBBBBBBBBBBBBBBBBBBBBBBBBBB");
			        }*/
			   //  Generamatrix=  vectorToMatrix( colonne1);
			    

			        double [][]SobPoints= Sob.getSobolPointSetForGeneratMatrix(d+1, w, Generamatrix );//ici k=1    d ou d+1???
			  
			       
			        
			        
			        
			        
			        Wafom lowWaf=new Wafom (Sob,1,dim,w,k);
			        
				       double currentWafom =lowWaf.calcWafom1(SobPoints,d+1);
				       
				       if (currentWafom < smallestWafom) {
				           smallestWafom = currentWafom;
				           bestColonne= colonne1; // Mise à jour de la colonne correspondante
				           //System.out.println(i1+" fois de plus");
				           
				           WAFOMS[d]=smallestWafom;
					    	bestColonnes[d]=bestColonne;
			 
				          // Generamatrix=replaceColumn(Generamatrix,bestColonne, d);
				       }

			        // Ajouter une nouvelle colonne à la matrice à chaque itération
			        //matrix = concatenateMatrixAndColumn(matrix, colonne1);
/*
			        // Affichage de la matrice résultante
			        System.out.println("Afficher Matrix");
			        for (int i = 0; i < matrix.length; i++) {
			            for (int j = 0; j < matrix[0].length; j++) {
			                System.out.print(matrix[i][j] + " ");
			            }
			            System.out.println();
			        }
			        
			        
			        */
			        
			        
			    	}
			    	
			    	
			    	}
			    	
			    	
			    	
			    	
			    	
			    	
			    	System.out.println("WAFOMS");
			    	System.out.println(Arrays.toString( WAFOMS));
			    	
			    	
			        System.out.println("BEST COLOMNS");
			        for (int i = 0; i <  bestColonnes[0].length; i++) {
			            for (int j = 0; j <  bestColonnes.length; j++) {
			                System.out.print(bestColonnes[j][i] + " ");
			            }
			            System.out.println();
			        }
			        
			    //}
			    
			    	
			        System.out.println("MAtrix");
			        for (int i = 0; i <  Generamatrix.length; i++) {
			            for (int j = 0; j <  Generamatrix[0].length; j++) {
			                System.out.print( Generamatrix[i][j] + " ");
			            }
			            System.out.println();
			        }
			    	
			 
	 
		  
			        
			        

			    	long endTime = System.currentTimeMillis();
			    	long executionTime = endTime - startTime;

			    	System.out.println("Temps d'exécution : " + executionTime + " millisecondes");
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

	 
}