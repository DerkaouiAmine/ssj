package fichier_Sob;

import java.io.FileNotFoundException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

public class TesterMeilleurWaf {
	
	 public static void main(String[] args) throws FileNotFoundException { 
		 long startTime = System.currentTimeMillis();


	int dim=3;
    int dimension=0;
	int k=14;//k=23 est le maximum dans notre cas
	int w=30;
	
	
	 double functionValue=0;
	 
	 
	 String chemin="/home/derkaoui/Desktop/ResultatsSob";
		RandomStream randomStream = new MRG32k3a();
		
		
		Tally stats=new Tally("Simul");
		Tally stats1=new Tally("Vari");
		
		int nbreSimulateLMS=100;
		int nbreSimulateDigi=100;
		//int N=(int) Math.pow(2, k);//equivalent a 1<<k
		
		double []moyenne = new double[nbreSimulateLMS] ;
		double [] WafomResult=new double[nbreSimulateLMS];
		
		int[][] []tableauPrincipalGenerLMS = new int[nbreSimulateLMS][][];
		
		double [][] pointSet=new double[1<<k][dim];
		
		//la plus bonne
		
	  /* int[][] Generamatrix = {
	            {1046994105, 494585195, 227311910, 113002895, 41890699, 26223387, 9376982, 5018071, 3480865, 1831082, 842385, 378132, 232991, 79423, 57829, 32287, 15145, 6217, 2313, 1270, 658},
	            {1069297368, 930578335, 619289281, 847556527, 708598175, 843628063, 828388775, 626005063, 821772494, 1037509775, 840876133, 772825372, 1008169578, 540187829, 688573138, 905776479, 572863156, 654057000, 922117121, 746174947, 671251987},
	            {1065712118, 536711880, 737016817, 556264442, 331072214, 939752460, 953392092, 396784615, 850592168, 597744856, 524525079, 551966362, 828597301, 355879173, 898039098, 696784456, 396152808, 593107583, 942327054, 352994066, 1017271835}
	        };*/
		
		//la plus mauvaise
	    
	/* int[][] Generamatrix =  {
	            {544048359, 272122765, 198640308, 67516457, 61301733, 30386411, 14160461, 7504560, 2823100, 1768949, 599595, 271779, 136438, 103240, 56454, 25970, 9994, 5463, 3001, 1282, 941},
	            {554621823, 545713435, 787716810, 766630911, 931721561, 629661659, 782011062, 952502329, 788402387, 790632794, 1071866377, 928967911, 627636103, 704332383, 911739782, 542114734, 993810251, 986890851, 918438229, 675014683, 579750507},
	            {537385419, 281480881, 874046701, 1062906097, 335997384, 587694912, 1060825998, 306951523, 1059965453, 1025534446, 352314123, 817873308, 709577199, 378552876, 850292664, 819922537, 270514576, 761689649, 779686648, 458286946, 966706506}
	        };*/
	 int[][] Generamatrix =  {
	 
	 {1063714292, 534518442, 263181785, 131568814, 62780033, 30924374, 13482947, 7839359, 3451979, 2050009, 859219, 518256, 201215, 95877},
			 {1068361547, 552232653, 1062257147, 765279082, 881589079, 980562997, 778556161, 664636942, 763600125, 587054550, 900105256, 591859850, 605455589, 1039282866},
			 {1069723551, 535901150, 537815274, 997758510, 345270997, 1018851399, 550925139, 504111212, 758938157, 895983491, 286876492, 671174418, 562284469, 471103001}
	 };
	 
	 
	    int vect[]=matrixToVector(Generamatrix);
	    
	    
	    
	    
	    
		for (int j=0;j<nbreSimulateDigi;j++) {
			stats.init();

		    DigitalNetBase2 Sob = new SobolSequence(Generamatrix[0].length,w,dim,vect);
			
			functionValue=0;
			Sob.addRandomShift(new MRG32k3a());
			
			
			pointSet=Sob.formatPointsTab();
				/*for (int a1=0;a1<dim;a1++)
				{
					double[] p=ExtractVecteur(pointSet,a1);
					//Methode 1
					functionValue*=evaluateFunction(p);
					//Methode2
					//functionValue*=(1 + randomStream.nextDouble() * (calculateNorm(p) - 0.5));
					
				}*/
			/*for (int i1 = 0; i1 < pointSet.length; i1++) {
			    double product = 1; // Déplacer l'initialisation ici, à l'intérieur de la boucle externe
			    RandomStream randomStream1 = new MRG32k3a();
			    for (int i2 = 0; i2 < pointSet[0].length; i2++) {
			        product*= (1+randomStream1.nextDouble()*(pointSet[i1][i2]-0.5) );
			    }
			    
			    stats.add(product);
			}*/
			for (int i1 = 0; i1 < pointSet.length; i1++) {
			    double somme = 0; // Déplacer l'initialisation ici, à l'intérieur de la boucle externe
			    RandomStream randomStream1 = new MRG32k3a();
			    for (int i2 = 0; i2 < pointSet[0].length; i2++) {
			        somme += Math.pow(randomStream1.nextDouble(), 2)* Math.pow((pointSet[i1][i2] - randomStream1.nextDouble()),2);
			    }
			    functionValue = Math.exp(-somme);
			   stats.add(functionValue);
			}
			
		    double variance=stats.variance();

			System.out.println(variance);
			
			
			
		/*	
			int [][]genera=new int[1][2];
			
			for(int i=0;i<genera.length;i++) {
				for(int j1=0;j1<genera[0].length;j1++) {
					
					System.out.print(genera[i][j1]+" ");
				}
				System.out.println();
			}*/
			
			
}
		
	
		
		
		
	 
	 
	 
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
