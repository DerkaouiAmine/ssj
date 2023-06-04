package fichier_Sob;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

public class SimulateSobol {
	 public static void main(String[] args) throws FileNotFoundException { 

	int dim=12;
    int dimension=0;
	int k=14;
	int w=31;
	
	
	 double functionValue=1;
	 
	String chemin="/home/derkaoui/Desktop/ResultatsSob";
	RandomStream randomStream = new MRG32k3a();
	
	
	Tally stats=new Tally("Simul");
	Tally stats1=new Tally("Vari");
	
	int nbreSimulateLMS=100;
	int nbreSimulateDigi=100;
	
	double []moyenne = new double[nbreSimulateLMS] ;
	double [] variance = new double[nbreSimulateLMS];
	
	int[][] []tableauPrincipalGenerLMS = new int[nbreSimulateLMS][][];
	double [][] pointSet=new double[2^k][dim];
	
	
	
	for (int i=0;i<nbreSimulateLMS;i++) {
		stats.init();
		pointSet=new double[2^k][dim];
		functionValue=1;
		
		DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
		Sob.leftMatrixScramble(new MRG32k3a());
		tableauPrincipalGenerLMS[i]=  Sob.getGeneratorMatrices(k);
		
		
			for (int j=0;j<nbreSimulateDigi;j++) {
				Sob.addRandomShift(new MRG32k3a());
				
				
				pointSet=Sob.formatPointsTab();
					for (int a1=0;a1<dim;a1++)
					{
						double[] p=ExtractVecteur(pointSet,a1);
						//Methode 1
						functionValue*=evaluateFunction(p);
						//Methode2
						//functionValue*=(1 + randomStream.nextDouble() * (calculateNorm(p) - 0.5));
						
					}
				stats.add(functionValue);

	}
	moyenne[i]=stats.average();
    variance[i]=stats.variance();
	}
	
	System.out.println("*************MOYENNE********************");
	for(int i=0;i<moyenne.length;i++)
	{System.out.println(moyenne[i]);}
	System.out.println("*************Variance********************");
	for(int i=0;i<variance.length;i++)
	{System.out.println(variance[i]);
	stats1.add(variance[i]);}
	System.out.println("*********************************");
	
	int[] result=findMinMax(variance);
	
	
	
	System.out.println("Best gene Matri LMS");
	
	
	for (int i=0;i<tableauPrincipalGenerLMS[0].length;i++)
	{
    	for (int j=0;j<tableauPrincipalGenerLMS[0][0].length;j++)
    	{
    		System.out.print(tableauPrincipalGenerLMS[result[0]][i][j]+" ");
    	}

		System.out.println();
	}
	
	System.out.println("Worst gene Matri LMS");
	
	
	for (int i=0;i<tableauPrincipalGenerLMS[0].length;i++)
	{
    	for (int j=0;j<tableauPrincipalGenerLMS[0][0].length;j++)
    	{
    		System.out.print(tableauPrincipalGenerLMS[result[1]][i][j]+" ");
    	}

		System.out.println();
	}
	
	
	stats1.report();	 }
	 
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
	    }
	
	
	 

