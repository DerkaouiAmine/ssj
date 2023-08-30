package fichier_Sob;

import java.util.Arrays;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;





//POur s=1

public class WAFOMPointSet {
    public static void main(String[] args) {
       /*// double[] decimalNumbers = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.11,0.12,0.13,0.14,0.15,0.16,0.17};
        int decimalPlaces = 3; // Nombre de décimales à prendre en compte =w
        int dim=3;
        
        
        //tres important ici j'ai pris w=k ie on a 2^k points 
        
       // int[] binaryArray = convertDecimalToBinary(decimalNumbers, decimalPlaces);
      int[] binaryArray = {1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1};

        //int[] binaryArray = {1,0,0,0,1,0,0,1,0,0,1,0,1,1,0,1,1,0,0,0,1,0,0,1,1,0,1,1,0,1,0,1,1,0,1,1,1,1,1,1,1,1};

       // int[] sum = WFP1(binaryArray, decimalPlaces,dim);
        //int [] binaryArray1= {0,0,1,1,0,1};
        int[] sum = WFPTest(binaryArray, decimalPlaces,dim);
        //int[] sum1 = WFP(binaryArray, decimalPlaces);
        
        
       // System.out.println(binaryArray.length);
      //  System.out.println("Le tableau  : " + Arrays.toString(binaryArray));
       // System.out.println("Le tableau sum : " + Arrays.toString(sum1));

        System.out.println("Le tableau sum : " + Arrays.toString(sum));
        int[][] tableauBidimensionnel = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        int[] tableauUnidimensionnel = transformerTableau(tableauBidimensionnel);

        // Affichage du tableau unidimensionnel
        for (int i = 0; i < tableauUnidimensionnel.length; i++) {
            System.out.print(tableauUnidimensionnel[i] + " ");
        }
    */
    	
    	int dim=3;
      
    	int k=16;
    	int w=31;
    	
    	
    	for (int i=0;i<100;i++) {
    	DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
		Sob.leftMatrixScramble(new MRG32k3a());
		Sob.addRandomShift(new MRG32k3a());
		//System.out.println(Sob.formatPoints());
		double[][] Tab= Sob.formatPointsTab();
//Sob.printGeneratorMatricesTrans(dim);		
		//Sob.printGeneratorMatrices(dim);
		//System.out.print(Sob.formatPoints());
		double []pointSet=transformerTableau(Tab);
		
		
	//	System.out.println(pointSet.length);
		
		/*for(int i=0;i<pointSet.length;i++)
		{
			System.out.print(pointSet[i]+"   ");
			
		}
		System.out.println();*/
		
		int[] binaryArray = convertDecimalToBinary(pointSet, w);
		//int[] binaryArray1= {1,0,0,0,0,1,0,1,0,1,1,1};
		//System.out.println(Arrays.toString(binaryArray));
		
		
		//System.out.println(Arrays.toString(binaryArray));
		//int[] binaryArray1= {1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1};
	
		
		 int[] sum=Wafom(binaryArray ,w,dim);
		 
		// System.out.println(Arrays.toString(sum));

		 //int[] sum=WFPTest2(binaryArray ,w,dim);
		
	   //   System.out.println("Le tableau  : " + Arrays.toString(binaryArray));

	     //  System.out.println("Le tableau sum : " + Arrays.toString(sum));
	       
	       //faire la somme
	       double tempsum=0;
	       double[] partialSum = new double[sum.length];
	       for (int r=0;r<sum.length;r++)
	       {
	    	   tempsum=tempsum+Math.pow(2, -sum[r]);
	    	   partialSum[r]= Math.pow(2, -sum[r]);
	       }
	       
	      // System.out.println(Arrays.toString(partialSum));
	       System.out.println(Math.log(tempsum));
    	}
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

    public static int[] WFP(int[] binaryArray, int decimalPlaces) {
        int[] sum = new int[binaryArray.length/decimalPlaces];
        int startIndex = 0;
        
        for (int i = 0; i < binaryArray.length/decimalPlaces; i++) {
            int tempSum = 0;
          
            for (int j = startIndex*decimalPlaces,k=1; j< ((startIndex+1)*decimalPlaces); j++,k++) {
 
                tempSum += k * binaryArray[j];
            }
            sum[i] = tempSum;
            startIndex++;
        }

        return sum;
    }

/*
public static int[] WFP2(int[] binaryArray, int decimalPlaces,int dim) {
    int[] sum = new int[binaryArray.length/decimalPlaces];
    int startIndex = 0;
    
    for (int i = 0; i < sum.length; i++) {
        int tempSum = 0;
      for (int l=0;l<dim;l++) {
    	  for (int j = startIndex*decimalPlaces,k=1; j< ((startIndex+1)*decimalPlaces); j++,k++) {
    		  
              tempSum += k * binaryArray[j];
          }
        startIndex++;
      }
        sum[i] = tempSum;
      
    }

    return sum;
}












public static int[] WFP1(int[] binaryArray, int decimalPlaces, int dim) {
    int[] sum = new int[binaryArray.length / decimalPlaces];
    int startIndex = 0;

    for (int i = 0; i < sum.length; i++) {
        int tempSum = 0;
        for (int l = 0; l < dim; l++) {
            for (int j = startIndex * decimalPlaces, k = 1; j <= (startIndex + 1) * decimalPlaces; j++, k++) {
                tempSum += k * binaryArray[j];
            }
           
            
        }
        sum[i] = tempSum;
        startIndex++;
        
    }

    return sum;
}






//prend par colonne
public static int[] WFPTest1(int[] binaryArray, int decimalPlaces, int dim) {
    int[] PartielSum = WFP(binaryArray, decimalPlaces);
    int index = PartielSum.length / dim,l=0;

   // System.out.println(PartielSum.length+ "  "+index);
    int[] sum = new int[index];
    int tempSum = 0;
   
	
    for (int i = 0; i < index; i++) {
    	 tempSum = 0;
		 l=0;
    	 for (int j = 0; j < dim; j++)
    	 {	
    		 tempSum += PartielSum[l + i];
    		 l+=index;
    		 
    	 }
        
       
        
        sum[i] = tempSum;
    }
    
    return sum;
}
*/


//prend par Ligne
public static int[] Wafom(int[] binaryArray, int decimalPlaces, int dim) {
  int[] PartielSum = WFP(binaryArray, decimalPlaces);
  int index = PartielSum.length / dim,l=0;

 // System.out.println(PartielSum.length+ "  "+index);
  int[] sum = new int[index];
  int tempSum = 0;
 
	
  for (int i = 0; i < index; i++) {
  	 tempSum = 0;
		 
  	 for (int j = 0; j < dim; j++)
  	 {	
  		 tempSum += PartielSum[j+l];
  		
  		 
  	 }
      l+=dim;
     
      
      sum[i] = tempSum;
  }
  
  return sum;
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

public static int produitScalaire(int[] vecteur1, int[] vecteur2) {
    if (vecteur1.length != vecteur2.length) {
        throw new IllegalArgumentException("Les vecteurs doivent avoir la même taille");
    }

    int produit = 0;
    for (int i = 0; i < vecteur1.length; i++) {
        produit += vecteur1[i] * vecteur2[i];
    }

    return produit;
}


}




