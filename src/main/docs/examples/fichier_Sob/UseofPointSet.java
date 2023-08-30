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
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;






//reste a teste le biais pour a variance  on prend encore lms +100 digital shift et encore 100 digital shift on vois si sa change ou pas pour voir si y'
public class UseofPointSet{
	 public static void main(String[] args) throws FileNotFoundException { 
		 long startTime = System.currentTimeMillis();


	int dim=4;
    int dimension=0;
	int k=3;
	int w=3;
	
	
	

		
		DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
		
	System.out.println(Sob.formatPoints());
double [][]tab=Sob.getSobolPointSet( k,w);
	


double [][]pointSet1=Sob.formatPointsTab();
int [][]pointSet=convertDecimalToBinary(pointSet1,w);

for (int i=0;i<tab.length;i++) {
	
	for (int j=0;j<tab[0].length;j++) {
		System.out.print(tab[i][j]+" ");
		
	}
	System.out.println();
	}


	 }
	 
	 
	 public static int[][] convertDecimalToBinary(double[][] decimalNumbers, int decimalPlaces) {
	        int[][] binaryArray = new int[decimalNumbers.length][decimalNumbers[0].length * decimalPlaces];

	        for (int i = 0; i < decimalNumbers.length; i++) {
	            for (int j = 0; j < decimalNumbers[i].length; j++) {
	                double decimalNumber = decimalNumbers[i][j];
	                int[] binaryRow = new int[decimalPlaces];

	                for (int k = 0; k < decimalPlaces; k++) {
	                    decimalNumber *= 2;
	                    binaryRow[k] = (int) decimalNumber;
	                    decimalNumber -= binaryRow[k];
	                }

	                System.arraycopy(binaryRow, 0, binaryArray[i], j * decimalPlaces, decimalPlaces);
	            }
	        }

	        return binaryArray;
	    }
	    }
	
	
	 

