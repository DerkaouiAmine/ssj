package fichier_Sob;

import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;

public class BROUILLONDUWAFOMS {
	 public static void main(String[] args) throws FileNotFoundException { 
		 int dim=3;
		    int dimension=0;
			int k=5;
			int w=5;
			//je vais le modifier pour agir directement sur le tableau principale
			int nombdeDeColonnesExtraites =1 ;
			
			
			
			int nbreSimulateLMS=100;
			int nbreSimulateDigi=100;
			
			double[][]WAFOMS=new double[nbreSimulateLMS][nbreSimulateLMS];
			int[][] []tableauPrincipalGenerLMS = new int[nbreSimulateLMS][][];
			
			int [] bestColonnes=new int [ nbreSimulateLMS];

		        
		        
		        
				
				//int [][] matrix = new int[dim][k];
				int[][] matrix = new int[dim][0]; // Matrice de sortie avec zéro colonne initialement

			    //for (int l = 0; l < k; l++) {//l=colonne extraite
			    	
			    	
				double smallestWafom=(1<<3);
			    	for(int i1=0;i1<100;i1++) {
			    	DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
					
					
					
					tableauPrincipalGenerLMS[0] =  Sob.getGeneratorMatrices(k);//ici sa sera un i au lieu de 0
			    	
			    	
			    	
			        int[] colonne1 = extraireColonne(tableauPrincipalGenerLMS[0], 0);//ici i1 a la place de 0
			        //matrix = concatenateMatrixAndColumn(matrix, colonne1);
			     //  double [][]tab= Sob.getSobolPointSetForGeneratMatrix(k, w, matrix);
			        Wafom lowWaf=new Wafom (Sob,1,dim,w,1);
				       double currentWafom =lowWaf.calcWafom();
				       
				       if (currentWafom < smallestWafom) {
				           smallestWafom = currentWafom;
				           bestColonnes= colonne1; // Mise à jour de la colonne correspondante
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
			    	System.out.println(Arrays.toString( bestColonnes)+" "+smallestWafom+"  ");
			        
			    //}
			    
			    
			    
			    
			    
			    
	 
		        
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
	 
}