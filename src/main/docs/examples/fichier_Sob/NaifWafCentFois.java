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
public class NaifWafCentFois {
	 public static void calc(int k) { 
		 long startTime = System.currentTimeMillis();
System.out.println("je calcule pour k="+k);
	int dim=3;
    int dimension=0;
	//int k=14;//k=23 est le maximum dans notre cas
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
	
	
	
	for (int i=0;i<nbreSimulateLMS;i++) {
		stats.init();
		pointSet=new double[1<<k][dim];
		
		
		DigitalNetBase2 Sob = new SobolSequence(k,w,dim);
		Sob.leftMatrixScramble(new MRG32k3a());
		tableauPrincipalGenerLMS[i] =  Sob.getGeneratorMatrices(k);
		Wafom waf=new Wafom (Sob,1,dim,w,k);
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
		 //   System.out.println("je suis a la "+i);
			
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
	int [] result1=findMinMax(WafomResult);
	double [] rearrangeVariance=arrangeVecteur(variance);
	double []rearrangeWafom = arrangeVecteur(WafomResult);
	
	System.out.println("Le plus petit Wafom : "+ WafomResult[result1[0]]);
	System.out.println("Le plus grand Wafom : "+ WafomResult[result1[1]]);
	
	System.out.println("Best gene Matri LMS");
	
	
	for (int i=0;i<tableauPrincipalGenerLMS[0].length;i++)
	{
    	for (int j=0;j<tableauPrincipalGenerLMS[0][0].length;j++)
    	{
    		System.out.print(tableauPrincipalGenerLMS[result1[0]][i][j]+" ");
    	}

    	System.out.println();
	}
	//System.out.println("bestWafom:"+WafomResult[result1[0]]);

	//System.out.println("LE rang de la meileur matrice: "+result[0]+"rang du meilleur wafom: "+result1[0]);
	
	System.out.println("Worst gene Matri LMS");
	
	
	for (int i=0;i<tableauPrincipalGenerLMS[0].length;i++)
	{
    	for (int j=0;j<tableauPrincipalGenerLMS[0][0].length;j++)
    	{
    		System.out.print(tableauPrincipalGenerLMS[result1[1]][i][j]+" ");
    	}

		System.out.println();
	
	}
	
	//System.out.println("worstWafom:"+WafomResult[result1[1]]);

	//System.out.println("LE rang de la worst matrice: "+result[1]+"rang du worst wafom: "+result1[1]);
	stats1.report();	
	
	
	 //System.out.println("Le meilleur wafom et son rang: "+WafomResult[result1[0]]+"----"+result1[0]);
	 //System.out.println("Le meilleur wafom selon le critere de variance: "+WafomResult[result[0]]+"----"+result[0]);
	 
	 //System.out.println("Le meilleur variance et son rang: "+variance[result[0]]+"----"+result[0]);
	 
	 //System.out.println("Le meilleur variance selon le critere wafom: "+variance[result1[0]]+"----"+result1[0]);
/*
System.out.println();
System.out.println("Le meilleur Wafom est : "+WafomResult[result1[0]]+" et a comme variance :"+variance[result1[0]]);
System.out.println("Le worst Wafom est : "+WafomResult[result1[1]]+" et a comme variance :"+variance[result1[1]]);
System.out.println("FAUT LE REVOIR CAR LES VARIANCE ");
System.out.println("Le wafom reduit la variance a l'ordre de "+(variance[result1[1]]/variance[result1[0]]));

double [] variance1=arrangeVecteur(variance);

System.out.println("En realite la plus petit est :"+variance1[0]+" le plus grand "+variance1[variance1.length-1]+" on a un ordre de :"+(variance1[variance1.length-1]/variance1[0]));
*/
	
//wafom du plus petit au plus grand
	
	int [] indicesTri=getIndicesTri(WafomResult);
	
	//System.out.println("Les 20 meilleur Wafom:");
	//je ppuvais utiliser juste les indices
	/*for(int i=0;i<20;i++) {
		//System.out.println(rearrangeWafom[i] +" variance correspanante : "+variance[indicesTri[i]]);
		System.out.println(WafomResult[indicesTri[i]]);

	}*/
	
	System.out.println("variance des  meilleur Wafom:");

	for(int i=0;i<(10);i++) {
		//System.out.println(rearrangeWafom[i] +" variance correspanante : "+variance[indicesTri[i]]);
		System.out.println(variance[indicesTri[i]]);

	}
	
	//System.out.println("Les 20 mauvaises Wafom:");
	System.out.println("variance des mauvaises Wafom:");

	for(int i=rearrangeWafom.length-1;(rearrangeWafom.length-(10))<=i;i--) {
		//System.out.println(rearrangeWafom[i] +" variance correspanante : "+variance[indicesTri[i]]);
		System.out.println(variance[indicesTri[i]]);

	}
	 


	long endTime = System.currentTimeMillis();
	long executionTime = endTime - startTime;

	System.out.println("Temps d'exécution : " + executionTime + " millisecondes");
	
 }
	 public static void main(String[] args) throws FileNotFoundException { 
		 
		 
		 for (int k=1;k<23;k++) {
			 calc(k);
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
	
	
	 

