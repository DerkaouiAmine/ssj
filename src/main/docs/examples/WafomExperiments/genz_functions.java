package WafomExperiments;

/**
 * In this class, we define our test functions as the functions of Genz, as defined in:
 * Alan Genz. “A Package for Testing Multiple Integration Subroutines”. In: 1987.
 */

//Oscillatory function: cos(2*pi*u1 + sum(a[i] * x[i]))
public class genz_functions {
	static double oscillatory(double[] x, double[] a,double [] u) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			sum += a[i] * x[i];
		}
		return Math.cos(2 * Math.PI *u[0]  + sum);
	}

	// Product Peak function: prod(a[i]^-2 + (x[i] - u[i])^-2)
	static double productPeak(double[] x, double[] a, double[] u) {
		double product = 1.0;
		for (int i = 0; i < x.length; i++) {
			double Result=1/(Math.pow(a[i], -2)+Math.pow(x[i]-u[i],2));

			product *= Result;
		}
		return product;
	}

	// Corner Peak function: (1 + sum(a[i] * x[i]))^(-s + 1)
	static double cornerPeak(double[] x, double[] a) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			sum += a[i] * x[i];
		}
		return Math.pow(1 + sum, -x.length + 1);
	}

	// Continuous function: exp(-sum(a[i] * |x[i] - u[i]|))
	static double continuous(double[] x, double[] a,double [] u) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			sum += a[i] * Math.abs(x[i] - u[i]);
		}
		return Math.exp(-sum);
	}

	// Discontinuous function: 0 if x[i] > u[i], otherwise exp(sum(a[i] * x[i]))
	static double discontinuous(double[] x, double [] a,double[] u) {
		for (int i = 0; i < x.length; i++) {
			if ((x[1] > u[1]) ||(x[2] > u[2])) {
				return 0.0;
			}
		}
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			sum += a[i] * x[i];
		}
		return Math.exp(sum);
	}


	static double gaussian(double[] x, double[] a, double[] u) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			sum += a[i] * a[i] * Math.pow(x[i] - u[i], 2);
		}
		return Math.exp(-sum);
	}
}
