package WafomExperiments;


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

public class WAfomExempleFinanceS16BIS {

	AsianOption asian;
	RandomStream noise = new MRG32k3a();



	public WAfomExempleFinanceS16BIS(AsianOption asian) {
		this.asian = asian;
	}



	// Main program: QMC and RQMC experiment with Asian option.
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
		asian.setProcess(new GeometricBrownianMotion(s0, r, sigma,
				new BrownianMotion(0, 0, 1, gen)));
		WAfomExempleFinanceS16BIS test = new WAfomExempleFinanceS16BIS(asian);

		Tally statValueMC = new Tally("Stats on payoff with crude MC");
		int n = 100000; // 10 million runs for Monte Carlo.
		Chrono timer = new Chrono();



		int s = dim;

		int w=30;


		for (int methode=0;methode<2;methode++) {
			System.out.println("******************Methode MCVSRQMC numero"+methode+"****************");



			WafomsStorage storage = new WafomsStorage(s,methode,1);
			int[][] Generamatrix = storage.getWafomsTabl();


			Tally stats1=new Tally("Simul RQMC");

			Tally statsrqmc=new Tally();


			int replicates=1000;



			for (int k1=10;k1<=Generamatrix[0].length;k1++) {
				// System.out.println("--------------------------pour k= "+k1);




				int [][] generat=extractSubMatrix(Generamatrix,Generamatrix.length-1,k1-1);


				int []vect=matrixToVector(generat);

				DigitalNetBase2 Sob = new SobolSequence(k1,w,dim,vect);
				stats1.init();

				PointSetRandomization rand = new RandomShift(new MRG32k3a());






				//  System.out.println(  RQMCExperiment. makeComparisonExperimentMCvsRQMC (asian, new MRG32k3a(), Sob, rand, 1000000, replicates)); 
				RQMCExperiment. simulReplicatesRQMC(asian, Sob, rand, replicates,stats1) ;

				System.out.println(stats1.average()+","+stats1.variance());


			}



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

	public static int[][] extractSubMatrix(int[][] originalMatrix,  int endRow , int endCol) {
		int startCol=0;
		int startRow=0;
		// Vérifier les limites des indices pour éviter les erreurs
		int numRows = originalMatrix.length;
		int numCols = originalMatrix[0].length;

		if (startRow < 0 || endRow >= numRows || startCol < 0 || endCol >= numCols) {
			throw new IllegalArgumentException("Les indices de sous-matrice sont hors limites.");
		}

		// Calculer la taille de la sous-matrice
		int subMatrixRows = endRow - startRow + 1;
		int subMatrixCols = endCol - startCol + 1;

		// Créer la sous-matrice
		int[][] subMatrix = new int[subMatrixRows][subMatrixCols];

		// Copier les éléments de la matrice originale dans la sous-matrice
		for (int i = startRow; i <= endRow; i++) {
			for (int j = startCol; j <= endCol; j++) {
				subMatrix[i - startRow][j - startCol] = originalMatrix[i][j];
			}
		}

		return subMatrix;
	}






}
