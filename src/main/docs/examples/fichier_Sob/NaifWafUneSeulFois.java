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
public class NaifWafUneSeulFois {
	 public static void main(String[] args) throws FileNotFoundException { 
		 long startTime = System.currentTimeMillis();


	int dim=3;
    int dimension=0;
	int k=23;//k=23 est le maximum dans notre cas
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
	double [] variance = new double[nbreSimulateLMS];
	double [] WafomResult=new double[nbreSimulateLMS];
	
	int[][] []tableauPrincipalGenerLMS = new int[nbreSimulateLMS][][];
	
	double [][] pointSet=new double[1<<k][dim];
	
	
	
	for (int i=1;i<=k;i++) {
		System.out.println("Je suis a la :"+i);
		
		stats.init();
		pointSet=new double[1<<k][dim];
		
		
		DigitalNetBase2 Sob = new SobolSequence(i,w,dim);
		Sob.leftMatrixScramble(new  MRG32k3a() );
		
		int [][]tab =  Sob.getGeneratorMatrices(i);
		for (int []v:tab)
		{
			System.out.println(Arrays.toString(v));
		}
		
		Wafom waf=new Wafom (Sob,1,dim,w,i);
		//System.out.println(waf.calcWafom());
		 WafomResult[i]=waf.calcWafom();
		
			for (int j=0;j<nbreSimulateDigi;j++) {
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
				
				
	}
			moyenne[i]=stats.average();
		    variance[i]=stats.variance();	
			//System.out.println("Le  Wafom : "+ WafomResult[i]);
			System.out.println(WafomResult[i]);
			
			
		

			

			

			
	}
	

	
}
	 
	 public static double[] ExtractVecteur(double[][] resultatMatr,int j)
		
		{double[] vj = new double[resultatMatr.length];
			for (int i = 0; i < resultatMatr.length; i++) {
	            vj[i] = resultatMatr[i][j];
	        }
			return vj;
		}
		
	 
	 
	 
		public static double evaluateFunction(double[] point) {
	        
	    	RandomStream randomStream = new MRG32k3a();
	    	double Somme=0;
	    	for (int i=0;i<point.length;i++)
	    	{
	    		Somme+=(1 + randomStream.nextDouble() * (point[i] - 0.5));
	    	}
	        return Somme/point.length;
		
		}
		
	    public static double calculateNorm(double[] vector) {
	        double sumOfSquares = 0.0;
	        
	        for (double component : vector) {
	            sumOfSquares += Math.pow(component, 2);
	        }
	        
	        return Math.sqrt(sumOfSquares);
	    }
	    
	    public static int[] findMinMax(double [] v)
	    {
	        int minIndex = 0;
	        int maxIndex = 0;
	        double minValue = v[0];
	        double maxValue = v[0];
	        
	        for (int i = 1; i < v.length; i++) {
	            if (v[i] < minValue) {
	                minValue = v[i];
	                minIndex = i;
	            }
	            
	            if (v[i] > maxValue) {
	                maxValue = v[i];
	                maxIndex = i;
	            }
	        }
	        int[] indices = {minIndex, maxIndex};
	        return indices;
	    }
	    
	    public static double[] arrangeVecteur(double[] vecteur) {
	        // Copie du vecteur original pour éviter de modifier l'original
	        double[] vecteurCopie = Arrays.copyOf(vecteur, vecteur.length);

	        // Tri du tableau accendant le 1er element est le plus petit le dernier et le plus grand
	        Arrays.sort(vecteurCopie);

	        return vecteurCopie;
	    }
	    

	    public static int[] getIndicesTri(double[] vecteurOrigi) {
	    	int []indices=new int[vecteurOrigi.length];
	    	 double [] vecteurOrdonne=arrangeVecteur(vecteurOrigi);
	    //	double []vecteur1=Arrays.copyOf(vecteur, vecteur.length);
	    	for(int i=0;i<vecteurOrigi.length;i++) {
	    		for(int j=0;j<vecteurOrigi.length;j++) {
	    		if (vecteurOrdonne [i]==vecteurOrigi[j])
	    		{
	    			indices[i]=j;
	    		}
	    		}
	    		
	    		
	    	}
			return indices;
	    
	    }
	    }
	
	
	 

