package WafomExperiments;


import java.io.FileNotFoundException;
import java.util.Arrays;

import umontreal.ssj.discrepancy.Wafom;
import umontreal.ssj.discrepancy.WafomAcceler;
import umontreal.ssj.discrepancy.WafomRQMC;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.NiedSequenceBase2;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;


/**
 * In this class, we measure the difference in computation times between our WAFOM implementation method and the WAFOM Accel method proposed by Shin Harase in "A search for extensible low-WAFOM point sets".
 * In: Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357.
 * 
 * Remark: Its performance can be enhanced with point sets other than Sobol.
 */

public class testTimeExecuteWafomAndtWafomAccel {

	public static void main(String[] args) throws FileNotFoundException { 

		int dim=3;


		int w=30;

		int q=3;



		/*	for (int i=1;i<=25;i++) {
    	DigitalNetBase2 Sob=new SobolSequence(i,w,dim);
		Sob.leftMatrixScramble(new MRG32k3a());




		 WafomAcceler wafAcce=	new  WafomAcceler (Sob, 1, dim,w, i, q);
		 Wafom waf=	new  Wafom (Sob, 1, dim,w, i);
	    	long startTime = System.currentTimeMillis();
	    	wafAcce.calcWafom();
	    	//wafAcce.calcWafom1(SobPoints, i);

	    	  	long endTime = System.currentTimeMillis();
	    	long executionTime = endTime - startTime;


	     	long startTime1 = System.currentTimeMillis();
	    	waf.calcWafom();
	    	//waf.calcWafom1(SobPoints,i);
	    	  	long endTime1 = System.currentTimeMillis();
	    	long executionTime1 = endTime1 - startTime1;

	    	System.out.println("****************temps d'excution pour k="+i+"**********************");

	    	System.out.println("Wafom Naif:"+ executionTime1 +" Wafom Accel:"+executionTime);

    	}*/


		// Avec Niedd

		for (int i=4;i<=16;i++) {

			DigitalNetBase2 Nied=new	NiedSequenceBase2 ( 20,  w,  i); 




			WafomAcceler wafAcce=	new  WafomAcceler (Nied, 1, i,w, 20, q);
			Wafom waf=	new  Wafom (Nied, 1, i,w, 20);
			long startTime = System.currentTimeMillis();
			wafAcce.calcWafom();

			long endTime = System.currentTimeMillis();
			long executionTime = endTime - startTime;


			long startTime1 = System.currentTimeMillis();
			waf.calcWafom();
			long endTime1 = System.currentTimeMillis();
			long executionTime1 = endTime1 - startTime1;


			System.out.println("****************temps d'excution pour d="+i+"**********************");

			System.out.println("Wafom Naif:"+ executionTime1 +" Wafom Accel:"+executionTime);

		}





	}


}