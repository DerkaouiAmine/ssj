package umontreal.ssj.discrepancy;

import java.io.FileNotFoundException;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;



//faut voir pour manipuler des blocs de byte au lieur de atbleau pour resuire la consomation de la memoire 
public class Wafom {
    private static final double factor = 1.0;
	private static DigitalNetBase2 dn;
	private static double c=1;//=par default
	private static int s;//dimension
	private static int w;
	private static int k;
    //private static final double factor = 2.0;
    public Wafom (DigitalNetBase2 dn, double c, int s,int w,int k) {
        this.dn=dn;
        this.c=c;
        this.s=s;
        this.w=w;
        this.k=k;
        
     }
    private static int m1p(int x) {
        return 1 - 2 * x;
    }

    private static double calcWafomSub(int[] point) {
        //int N = Long.SIZE;
    	int N=w;
        
        double prod = 1.0;
        int startIndex=0;

        for (int i = 0; i < s; i++) {
            //double w1 = point[i];
        	int index=startIndex*w;
            for (int j = 0; j < N; j++) {
                //int bij = (int) (((long)w >> (N - j - 1)) & 1);
            	int bij=point[j+index];
            	//System.out.println("Index:"+(j+index)+":BIJ:"+bij);
                double p = 1.0 + m1p(bij) * Math.pow(2.0, -factor * (j + 2 - c));
                prod *= p;
            }
            startIndex++;
            //System.out.println("AAAAAAAAAAAAAAAAAAAAAAAAA");
        }

        return prod - 1.0;
    }

    public static double calcWafom(){
        double sum = 0.0;
      
        
       // double [] partialProd=new double [(int) num];
        double[][] tab= dn.formatPointsTab();
       // double[][] tab= dn.getSobolPointSet(k, w);

		int [][]pointSet=convertDecimalToBinary(tab,w);
		long num = pointSet.length;
        for (long i = 0; i < num; i++) {
            int[] point = pointSet[(int) i];
          
            double sub = calcWafomSub(point);
            //partialProd[(int) i]=sub;
            sum += sub;
        }
//j'ai le produit partiel et je stocke a la fin la somme
        //partialProd[num]=sum/num;
        return sum/ num;
    }
    //Wafom partielle
    
    public static double calcWafom1(double[][] array,int k1) {
        double sum = 0.0;
        int num1=(1<<k1);
        int[][] pointSet = convertDecimalToBinary(array, w);
        long num = pointSet.length;
        for (long i = 0; i < num; i++) {
            int[] point = pointSet[(int) i];
            double sub = calcWafomSub(point);
            sum += sub;
        }
        return sum / num1;
    }

    
    
    public static int[][] convertDecimalToBinary(double[][] decimalNumbers, int decimalPlaces) {
        int[][] binaryArray = new int[decimalNumbers.length][decimalNumbers[0].length * decimalPlaces];

        for (int i = 0; i < decimalNumbers.length; i++) {
            for (int j = 0; j < decimalNumbers[i].length; j++) {
                double decimalNumber = decimalNumbers[i][j];
                int[] binaryRow = new int[decimalPlaces];

                for (int k = 0; k < decimalPlaces; k++) {
                    decimalNumber *= 2;
                    binaryRow[k] = (int) decimalNumber;
                    decimalNumber -= binaryRow[k];
                }

                System.arraycopy(binaryRow, 0, binaryArray[i], j * decimalPlaces, decimalPlaces);
            }
        }

        return binaryArray;
    }

  /*  public static void main(String[] args) throws FileNotFoundException { 
    	
    	int dim=5;
        
    	int k=16;
    	int w=31;
    	
    	
    	//for (int i=1;i<=25;i++) {
    	DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
		Sob.leftMatrixScramble(new MRG32k3a());
		//Sob.addRandomShift(new MRG32k3a());
		//System.out.println(Sob.formatPoints());
		double[][] tab= Sob.formatPointsTab();
		int [][]pointSet=convertDecimalToBinary(tab,w);
		
		/*for(int r1=0;r1<pointSet.length;r1++) {
		for(int j=0;j<pointSet[0].length;j++)
			{
				System.out.print(pointSet[r1][j]+"  ");
			}
			System.out.println("AAAAAAAAAA");
		}*/
    	//double a=calcWafom(pointSet, 1,w);
    	//System.out.println(Math.log10(a));
    	//}
   // }

}