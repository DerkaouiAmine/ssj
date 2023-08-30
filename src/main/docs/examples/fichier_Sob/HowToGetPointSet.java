package fichier_Sob;

import java.io.BufferedWriter;



import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import umontreal.ssj.hups.DigitalNet;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.stat.Tally;




// Faut l'ajouter a digitalNet




public class HowToGetPointSet{   
	//j =dim soit 0 ou 1 ou 2,..<s=dim
/*	public static double[][] calculateU(int s,int n, int k, int[][][] C1, int r) {
		// int n=(1 << k);
		
    double[][] U = new double[s][n];
    double normFactor = 1.0 / (1 << r);
   for(int j=0;j<s;j++) {
	   int [][]C=C1[j];
    for (int i = 0; i < n-1; i++) {
        int coord = 0; // Initialize coord as an integer
        for (int c = 0; c < k-1; c++) {
            coord ^= ((i >> c) & 1) * C[j][c]; // Convert C[j][c] to integer
            
        }
        U[j][i] = coord * normFactor;
    }
   }

    return U;
}*/
	
	public static double calculateU1(int j,int s,int n, int k, int[][] C, int r) {
		// int n=(1 << k);

		 double U = 0;
    double normFactor = 1.0 / (1 << r);
   
	   
    for (int i = 0; i < n; i++) {
        int coord = 0; // Initialize coord as an integer
        for (int c = 0; c < k; c++) {
           // System.out.println(c +" CCC");

            coord ^= ((i >> c) & 1) * C[j][c]; // Convert C[j][c] to integer
           //System.out.println(coord +" COOOR");
          
          // System.out.println(U[j][i]+" UUUUU");
        }
        
        U = coord * normFactor;
    }
   

    return U;
}
	
	

public static void main(String[] args) {
	/*int dim=3;
    
	int k=2;
	int w=5;
	
	
	
	DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
	System.out.print(Sob.formatPoints());;
	int n=(1<<k);
	
	int [][][]tab=Sob.generatorMatricesToStandardFormat();
	int [][]tab1=tab[0];
	for (int i=0;i<tab1.length;i++) {
		for (int j=0;j<tab1[0].length;j++)
		{System.out.print(tab1[i][j]+" ");}
		System.out.println();
	}
	

	
	double [][] a=calculateU(dim,n,  k, tab, w);
	
	for (int i=0;i<a.length;i++) {
		for (int j=0;j<a[0].length;j++)
		{System.out.print(a[i][j]+" ");}
		System.out.println();
	}
	/*double[][] U =new double[n][dim];
	
	int [][][]tab=Sob.generatorMatricesToStandardFormat();
	for (int i=0;i<tab.length;i++) {
		int [][]tab1=tab[i];
		for (int j=0;j<tab[0].length;j++) {
			for (int l=0;l<tab[0][0].length;l++) {
				U[n][i] = calculateU(dim,n, k, tab1, w);
		}
		
	}
		
	}
	
*/
	// System.out.println(Arrays.toString(U));*/

	int dim=3;
    
	int k=3;
	int w=7;
	
	
	
	DigitalNetBase2 Sob=new SobolSequence(k,w,dim);
	Sob.leftMatrixScramble(new MRG32k3a());;
	System.out.print(Sob.formatPoints());
	int n=(1<<k);
	
	int [][]tab=Sob.getGeneratorMatrices(k);
	

	

	double[][] a=new double[dim][n+1];
	for (int j=0;j<dim;j++) {
for (int i=1;i<n+1;i++) {
	
	a[j][i]=calculateU1(j,dim,i,  k, tab, w);
}
	}
	for (int j=0;j<dim;j++) {
for (int i=0;i<n+1;i++) {
	System.out.print(a[j][i]+" ");
}
	System.out.println();}
	

	
	
	
	
	
	
	int [] tab1=Sob.getGeneratorMatricesTrans();
			
		//	System.out.println(Arrays.toString(tab1));
			
			
			
			
			
			
}





}
