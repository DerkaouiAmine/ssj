package fichier_Sob;

import java.io.FileNotFoundException;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;

public class WafomTestAvecGener {
	 public static void main(String[] args) throws FileNotFoundException { 

		 int dim=3;
		    int dimension=0;
			int k=7;
			int w=7;
    
DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
		
Sob.leftMatrixScramble(new MRG32k3a());

		
		int [][] generat=  Sob.getGeneratorMatrices(k);//ici sa sera un i au lieu de 0
    	
    	
    	System.out.print(Sob.formatPoints());
        //matrix = concatenateMatrixAndColumn(matrix, colonne1);
    	
    	   System.out.println("Afficher genera");
	        for (int i = 0; i < generat.length; i++) {
	            for (int j = 0; j < generat[0].length; j++) {
	                System.out.print(generat[i][j] + " ");
	            }
	            System.out.println();
	        }
    	
	       int[][] tab2= Sob.getGeneratorMatrices(k);
    	
    	
      double [][]tab= Sob.getSobolPointSetForGeneratMatrix(2, w, tab2);//k=1,2 commence toujours par 1 generatr w*k y'aura au moins 1 colonne
      // Affichage de la matrice résultante
        System.out.println("Afficher Sobol");
        for (int i = 0; i < tab.length; i++) {
            for (int j = 0; j < tab[0].length; j++) {
                System.out.print(tab[i][j] + " ");
            }
            System.out.println();
        }
    
        Wafom lowWaf=new Wafom (Sob,1,dim,w,k);
	       double currentWafom =lowWaf.calcWafom1(tab);
	       System.out.println(currentWafom);
	       
	       double waf=lowWaf.calcWafom();
	       System.out.println(waf);
	       
	       
double a=0.03270964787696021;
double b=0.031383261689337624;

if (a>b) {System.out.print("VRAI");}
else {System.out.print("FAUX");}
	      

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

}
