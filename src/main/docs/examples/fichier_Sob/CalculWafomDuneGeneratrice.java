package fichier_Sob;

import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;

public class CalculWafomDuneGeneratrice {
	 public static void main(String[] args) {
		 
		 int dim=3;
			int k=14;
			int w=30;
		 
			/*int[][] Generamatrix =  {
					
					{1072724468, 536842282, 268372051, 134183710, 54487121, 24969208, 12400783, 7665377, 2676686, 2028212, 850310, 294855, 161447, 122743},
							{1073507273, 736034528, 713838568, 629559845, 551005190, 599263309, 1000066540, 944022191, 607990017, 600339377, 1031718804, 677633979, 823937335, 863404580},
							{1073297697, 536870813, 617464064, 768422308, 438458659, 640423542, 1006017563, 436752814, 816265633, 816182509, 463231848, 795433166, 740063652, 403703910}
				 };*/
			
			 int[][] Generamatrix =  {
					 {1073675440, 536720501, 268433701, 134215370, 66977650, 33551701, 16644777, 8385746, 4060028, 1965761, 1045787, 389326, 258126, 75486},
					 {1073701870, 537640631, 1073016086, 776690365, 906704984, 998279003, 758495926, 650189617, 898530890, 971159730, 637281622, 813755065, 985821855, 842834627},
					 {1073654377, 536836207, 536883239, 1000855111, 341533983, 1012493228, 542727479, 505936431, 883352971, 999362912, 282347294, 670829383, 556296526, 483125218}		
			 
			 };
					 
			int []vect=matrixToVector(Generamatrix);
			
			DigitalNetBase2 Sob = new SobolSequence(k,w,dim,vect);
			
			
			 Wafom lowWaf=new Wafom (Sob,1,dim,w,k);
		        
		       double currentWafom =lowWaf.calcWafom();
		       
		       System.out.println(currentWafom);
		       
		       
		       
		       
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
	 
	 
	  public static void growingVector1(int[] inputArray, int element, int place) {
	        int[] Vector = new int[inputArray.length];
	        inputArray[place]=element;
	        
	    }
	 
	 
    public static int[] growingVector(int initialElement, int iterations) {
        int[] growingVector = new int[iterations];

        // Add the initial element to the vector
        growingVector[0] = initialElement;

        // Iteratively add one element at each iteration
        for (int i = 1; i < iterations; i++) {
            int newElement = growingVector[i - 1] + 1;
            growingVector[i] = newElement;
        }

        return growingVector;
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
