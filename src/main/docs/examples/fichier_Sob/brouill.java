package fichier_Sob;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.stat.Tally;

public class brouill {
	 public static void main(String[] args) throws FileNotFoundException {
		 

			int dim=3;
		    int dimension=0;
			int k=14;
			int w=31;
			
			
			 double functionValue;
			String chemin="/home/derkaoui/Desktop/ResultatsSob";
			//RandomStream randomStream = new MRG32k3a();
			
			
			Tally stats=new Tally("Simul");
			
			
			int nbreSimulateLMS=1000;
			int nbreSimulateDigi=1000;
			
			double []moyenne = new double[nbreSimulateLMS] ;
			double [] variance = new double[nbreSimulateLMS];
			int[][][] tableauPrincipalGenerLMS = null;
			
			
			
			for (int i=0;i<nbreSimulateLMS;i++) {
				stats.init();
				
				DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
				Sob.leftMatrixScramble(new MRG32k3a());
				
				int[][] matrice=Sob.getGeneratorMatrices(k);
				
				tableauPrincipalGenerLMS[i]=matrice;
	
			//moyenne[i]=stats.average();
		    //variance[i]=stats.variance();
			}
			
			
			 System.out.print(tableauPrincipalGenerLMS[0]);
			 
			 
	 }
	 
	}
	
