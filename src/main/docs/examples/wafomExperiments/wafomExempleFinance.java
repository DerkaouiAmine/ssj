package wafomExperiments;

/**
 * In this class, we test our best point sets obtained in terms of WAFOM with an Asian option model.
 * We draw inspiration from the code defined in SSJ IFt6561.
 *  */

import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import ift6561examples.AsianOption;
import umontreal.ssj.hups.*;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.*;
import umontreal.ssj.mcqmctools.*;

public class wafomExempleFinance {

	AsianOption asian;
	RandomStream noise = new MRG32k3a();

	public wafomExempleFinance(AsianOption asian) {
		this.asian = asian;
	}

	public static void main(String[] args) {

		int dim = 16;

		int numObsTimes = dim;
		double T1 = 1.0 / (double) numObsTimes;
		double T = 1.0;
		double strike = 101.0;
		double s0 = 100.0;
		double r = 0.1;
		double sigma = 0.12136;
		RandomStream noise = new MRG32k3a();
		NormalGen gen = new NormalGen(noise);
		AsianOption asian = new AsianOption(r, numObsTimes, T1, T, strike);
		asian.setProcess(new GeometricBrownianMotion(s0, r, sigma, new BrownianMotion(0, 0, 1, gen)));
		wafomExempleFinance test = new wafomExempleFinance(asian);

		Tally statValueMC = new Tally("Stats on payoff with crude MC");
		int n = 100000; // 10 million runs for Monte Carlo.
		Chrono timer = new Chrono();

		int s = dim;

		int w = 30;

		for (int methode = 0; methode < 3; methode++) {
			System.out.println("******************Methode MCVSRQMC numero" + methode + "****************");
			int BestOrWorstWaf = 0;

			// we take only the best with this if
			/*
			 * if (methode==0||methode==1|| methode==2 ) { BestOrWorstWaf++; }
			 */

			while (BestOrWorstWaf < 2) {
				if (BestOrWorstWaf == 0) {
					System.out.println("******************Worst****************");
				}
				if (BestOrWorstWaf == 1) {
					System.out.println("******************Best****************");
				}

				WafomsStorage storage = new WafomsStorage(s, methode, BestOrWorstWaf);
				int[][] Generamatrix = storage.getWafomsTabl();

				Tally stats1 = new Tally("Simul RQMC");

				Tally statsrqmc = new Tally();

				int replicates = 1000;

				for (int k1 = 10; k1 <= Generamatrix.length; k1++) {
					// System.out.println("--------------------------pour k= "+k1);

					DigitalNetBase2 Sob = new DigitalNetBase2(k1, w, dim, Generamatrix[k1 - 1]);
					stats1.init();

					PointSetRandomization rand = new RandomShift(new MRG32k3a());

//					System.out.println(  RQMCExperiment. makeComparisonExperimentMCvsRQMC (asian, new MRG32k3a(), Sob, rand, 1000000, replicates)); 

					RQMCExperiment.simulReplicatesRQMC(asian, Sob, rand, replicates, stats1);

					System.out.println(stats1.average() + "," + stats1.variance());

				}

				BestOrWorstWaf++;
			}

		}

	}

}
