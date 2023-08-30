package fichier_Sob;
/*
 * L"idede est de mettre ;le tt dns un tale [] et le partager apres en snbre de dimension ie de 0 a 2^k puis de 2**K +1 juska 2*2^k...
 * */

import java.util.Arrays;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;

public class WAFOMbrouilloProduit { public static void main(String[] args) {
	int dim=3;
    
	int k=16;
	int w=31;
	
	
	for(int i=0;i<100;i++) {
	DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
	//Sob.leftMatrixScramble(new MRG32k3a());
	//Sob.addRandomShift(new MRG32k3a());
	//System.out.println(Sob.formatPoints());
	double[][] Tab= Sob.formatPointsTab();

	double []pointSet=transformerTableau(Tab);
	int[] binaryArray = convertDecimalToBinary(pointSet, w);
//	System.out.println(Arrays.toString(binaryArray));
//	int[] binaryArray1= {1,0,0,0,1,1};
	//int[] binaryArray1= {1,0,0,0,0,1,0,1,0,1,1,1};
	double wafom=calculateWAFOM(binaryArray,dim,w);
    System.out.println(Math.log(wafom));
	}
	
}

//idee tres importante faire un test pour voir si a^2-b^2 ou (a+-b)^2)


public static double calculateWAFOM(int[] pointSet, int dim, int w) {
    double norm = (pointSet.length/w)/dim;
    int factor=1;
    int c=1;
    double sum = 0;
    int startIndex = 0;
    
    for (int iter = 0; iter < norm; iter++) {
        double produit = 1.0;
        
        for (int j = 1; j <= dim; j++) {
            int index = startIndex * w;
          //  System.out.println("***********AAAAAAA"+index);
            for (int l = 1; l <= w; l++) {
                int a_ij = pointSet[index + l - 1];
               // int bij = (int) ((a_ij >> (w - j - 1)) & 1);
               // System.out.println(a_ij+"        "+bij);

                //System.out.println("***********"+(index+l-1));

                produit *= 1+(Math.pow(-1, a_ij)* Math.pow(2, -factor*(l+2-c)));
               
            }
            
            startIndex++;
        }
        
        sum +=(produit-1) ;
    }
    
    return  (sum/norm );
}


public static double calculateWAFOM1(int[] pointSet, int dim, int w) {
    double norm = (pointSet.length/w)/dim;
    
    double sum = 0;
    int startIndex = 0;
    
    for (int iter = 0; iter < norm; iter++) {
        double produit = 1.0;
        
        for (int j = 1; j <= dim; j++) {
            int index = startIndex * w;
            //System.out.println("***********AAAAAAA"+index);
            for (int l = 1; l <= w; l++) {
                int a_ij = pointSet[index + l - 1];
                int bij = (int) ((a_ij >> (w - j - 1)) & 1);
                //System.out.println("***********"+(index+l-1));

                produit *= 1+(Math.pow(-1, bij)* Math.pow(2, -l));
               
            }
            
            startIndex++;
        }
        
        sum +=(produit-1) ;
    }
    
    return  (sum/norm );
}

public static double[] transformerTableau(double[][] tableauBidimensionnel) {
    int nbLignes = tableauBidimensionnel.length;
    int nbColonnes = tableauBidimensionnel[0].length;
    int taille = nbLignes * nbColonnes;

    double[] tableauUnidimensionnel = new double[taille];
    int index = 0;

    for (int i = 0; i < nbLignes; i++) {
        for (int j = 0; j < nbColonnes; j++) {
            tableauUnidimensionnel[index] = tableauBidimensionnel[i][j];
            index++;
        }
    }

    return tableauUnidimensionnel;
}


public static int[] convertDecimalToBinary(double[] decimalNumbers, int decimalPlaces) {
    int[] binaryArray = new int[decimalNumbers.length * decimalPlaces];

    int index = 0;
    for (double decimalNumber : decimalNumbers) {
        // Conversion du nombre décimal en binaire
        double fractionalPart = decimalNumber - (int) decimalNumber;
        for (int i = 0; i < decimalPlaces; i++) {
            fractionalPart *= 2;
            binaryArray[index++] = (int) fractionalPart;
            fractionalPart -= binaryArray[index - 1];
        }
    }

    return binaryArray;
}
}
