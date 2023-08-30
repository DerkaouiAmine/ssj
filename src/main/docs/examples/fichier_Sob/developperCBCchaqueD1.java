package fichier_Sob;



import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;


//Methode2 ilfaut calculer le Waf pour chaque vecteur colonne et pas l'inserre a chaque fois pour avoir une matrice et comparerr les deux methodes

public class developperCBCchaqueD1 {//marche tres mal ie meme si on ades wafom tres perit pour des 1D sa veut pas dire que la matrice aura un petit Wafom
	 public static void main(String[] args) throws FileNotFoundException { 
	    	long startTime = System.currentTimeMillis();

		 int dim=3;
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
			
			int[] bestColonnes=new int [ k];
			int[] worstColonnes=new int [ k];


		        
		        
		        
				
				//int [][] matrix = new int[dim][k];
				int[]Generamatrix = new int[k]; // Matrice de sortie avec zéro colonne initialement

			    //for (int l = 0; l < k; l++) {//l=colonne extraite
			    	for(int dimension=1;dimension<=dim;dimension++) {
			    	for(int d=0;d<k;d++) {
				double smallestWafom=(1<<k);
				double biggestWafom=-1;
				//int[] colonne1=new int[dim];
			    	for(int i1=0;i1<1000;i1++) {
			    		
			    		
			    	DigitalNetBase2 Sob = new SobolSequence(d+1,w,dimension);
				       Sob.leftMatrixScramble(new MRG32k3a());
					
				       int [][] generatriOriginale=Sob.getGeneratorMatrices(d+1);
				       int element=generatriOriginale[dimension-1][d];
				       growingVector(Generamatrix, element, d);
				       int [] Generamatrix1 = Arrays.copyOfRange(Generamatrix, 0, d+1);
				       
				       DigitalNetBase2 Sob1 = new SobolSequence(d+1,w,1,Generamatrix1);//car je fais dimension par dimension donc 1 seul ce qui importe c'est les elemeents qui constitu la matrice generatrice

					      double [][]SobPoints = Sob1.formatPointsTab();
				     
					      Wafom lowWaf=new Wafom (Sob1,1,1,w,d+1);

					       double currentWafom =lowWaf.calcWafom1(SobPoints,d+1);
					       
					       if (currentWafom < smallestWafom) {
					    	   smallestWafom = currentWafom;
	    //System.out.println(i1+" fois de plus");
					    	   Generamatrix[d]=element;
					    	   bestWAFOMS[d]=smallestWafom;
					    	   	bestColonnes[d]=element;
					    	   	//il faut d*index;

					       }

					       if (currentWafom > biggestWafom) {
					    	   biggestWafom= currentWafom;
	   // bestColonne= colonne1; // Mise à jour de la colonne correspondante
	    //System.out.println(i1+" fois de plus");
	    
					    	   worstWAFOMS[d]=biggestWafom;
					    	   worstColonnes[d]=element;

					       }
				       

			    	}
			    	System.out.println(smallestWafom);
			    	
			    	
			    	}
			    	
			    	System.out.println(Arrays.toString(bestColonnes));

			    	
	 }
			  
	 
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