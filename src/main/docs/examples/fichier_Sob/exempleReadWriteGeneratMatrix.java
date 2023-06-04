package fichier_Sob;

import java.io.FileNotFoundException;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;

public class exempleReadWriteGeneratMatrix {
	 public static void main(String[] args) throws FileNotFoundException { 
	    	int dim=3;
	        int k=14;
	    	int w=31;
	    	
	    	String chemin="/home/derkaoui/Desktop/ResultatsSob";
	    	
	    	/*for (int i = 0; i < 100; i++) {
	    		Sob.leftMatrixScramble(new MRG32k3a());
	    	    int[][] tab = Sob.getGeneratorMatrices(k);
	    	    writeGeneratorMatricesToFile(chemin, tab);
	    	}*/


	    	int nbreSimulateLMS=1000;
	    	DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
	    	 int[][][] tableauPrincipalGenerLMS= new int[nbreSimulateLMS][][];
	    	//int [][] matrice=Sob.getGeneratorMatrices(k);
	    	//Sob.printGeneratorMatricesTrans(dim);
	    	
	    	

	    	for (int i = 0; i < nbreSimulateLMS; i++) {
	    		Sob.leftMatrixScramble(new MRG32k3a());
	    		//methode 1
	    	   Sob.writeGeneratorMatrixToFile(chemin);
	    	    
	    	  //GeneMatrix=Sob.readGeneratorMatrixFromFile(chemin);
	    	    
	    	    
	    	    
	    	    
	    	   
				//methode 2
	    	    //tableauPrincipalGenerLMS[i]=  Sob.getGeneratorMatrices(k);
	    	    
	    	    
	    	}

	    	

	    	
	    	
	    	 }
	    	
	    
	    

	    	
	    }

	
	



