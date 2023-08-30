package fichier_Sob;



import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;


//Methode2 ilfaut calculer le Waf pour chaque vecteur colonne et pas l'inserre a chaque fois pour avoir une matrice et comparerr les deux methodes

public class DevelopperCBCMArchepourD1 {
	 public static void main(String[] args) throws FileNotFoundException { 
	    	long startTime = System.currentTimeMillis();

		 int dim=1;
		    int dimension=0;
			int k=14;
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
			int [] bestColonne=new int[dim*k];
			
			int[] bestColonnes=new int [ k*dim];
			int[] worstColonnes=new int [ k*dim];


		        
		        
		        
				
				//int [][] matrix = new int[dim][k];
				int[]Generamatrix = new int[dim*k]; // Matrice de sortie avec zéro colonne initialement

			    //for (int l = 0; l < k; l++) {//l=colonne extraite
			    	int index=0;
			    	for(int d=index*k;d<((index+1)*k);d++) {
			    		index++;
				double smallestWafom=(1<<k);
				double biggestWafom=-1;
				//int[] colonne1=new int[dim];
			    	for(int i1=0;i1<100;i1++) {
			    		
			    		
			    	DigitalNetBase2 Sob = new SobolSequence(d+1,w,dim);
				       Sob.leftMatrixScramble(new MRG32k3a());
					
				       int [] generatriOriginale=Sob.getGeneratorMatricesTrans();

				       int element=extractElement(generatriOriginale,d);
				       growingVector(Generamatrix, element, d);

				       int [] Generamatrix1 = Arrays.copyOfRange(Generamatrix, 0, d+1);
				       


				       DigitalNetBase2 Sob1 = new SobolSequence(d+1,w,dim,Generamatrix1);

				      double [][]SobPoints = Sob1.formatPointsTab();

				       Wafom lowWaf=new Wafom (Sob1,1,dim,w,d+1);

				       double currentWafom =lowWaf.calcWafom1(SobPoints,d+1);

				       if (currentWafom < smallestWafom) {
				    	   smallestWafom = currentWafom;
    //System.out.println(i1+" fois de plus");
				    	   Generamatrix[d]=element;
				    	   bestWAFOMS[d]=smallestWafom;
				    	   	bestColonnes[d]=element;

				       }

				       if (currentWafom > biggestWafom) {
				    	   biggestWafom= currentWafom;
   // bestColonne= colonne1; // Mise à jour de la colonne correspondante
    //System.out.println(i1+" fois de plus");
    
				    	   worstWAFOMS[d]=biggestWafom;
				    	   worstColonnes[d]=element;

				       }

			    	}
		
			    	
			    	
			    	
			    	System.out.println("BEST COLOMNS");
			        System.out.println(Arrays.toString(bestColonnes));}
			    	long endTime = System.currentTimeMillis();
			    	long executionTime = endTime - startTime;
			    	
			    	//System.out.println(executionTime);
	 
			        

	 
			  
	 
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
}