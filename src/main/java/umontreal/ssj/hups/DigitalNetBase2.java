/*
 * Class:        DigitalNetBase2
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.hups;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.stream.IntStream;

import umontreal.ssj.rng.*;
import umontreal.ssj.util.*;

/**
 * A special case of @ref DigitalNet for the base @f$b=2@f$. The implementation exploit the binary
 * nature of computers and is much more efficient than for the case of a general @f$b=2@f$. Binary
 * expansions are easy to obtain because the computer already uses them internally. The generator
 * matrices @f$\mathbf{C}_j@f$ are stored in a single large array of size @f$sk@f$. The @f$c@f$-th column
 * of @f$\mathbf{C}_j@f$, for @f$c=0,\dots,k-1@f$, is stored at position @f$jk + c@f$ of this array,
 * as a 32-bit integer. For all derived classes, this 32-bit integer must have a binary 
 * representation of the form @f$
 * [0\, 0\, \cdots\, c_0\, c_1 \cdots \, c_{r-1}]@f$ (with the least significant bits on the right),
 * i.e., should be less than @f$2^r@f$. 
 * The value of @f$k@f$ cannot exceed 31 (32 is
 * not allowed because Java does not have 32-bit unsigned integers). 
 * The value of @f$w@f$ is always 31.
 * In this class, the random digital shift in base 2 generated by `addRandomShift` corresponds
 * to a random XOR.
 * 
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DigitalNetBase2 extends DigitalNet {
	private transient int[] originalMat;    // Original matrices, without randomization.
	protected int[] genMat;       // The current generator matrix.
	protected transient int[] digitalShift;   // Stores the digital shift vector.
	 //  **Pierre:** Check if @f$w@f$ is used sometimes in place of @f$r@f$, and clarify.
	private int[][] scrambleMat;
	private boolean fichierExiste = false;
	
	
	
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DERKAOUI ADD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/**
	 * Returns the generator matrices in decimal format, represented as a 2-dimensional array of shape dim x k.
	 * This format uses an int for each column of w bits, resulting in a more memory-efficient and faster method.
	 */

    public int[][] getGeneratorMatrices(int k) {
    	int[] GeneratorMatricesTrans= Arrays.copyOf(genMat, genMat.length);
    	int [][]GeneratorMatrices=new int[dim][k];
    	int l=0;
    	for (int i=0;i<dim;i++)
    	{l=0;
    	for(int j=i*k;j<(i+1)*k;j++)
    	{
    		GeneratorMatrices[i][l]=GeneratorMatricesTrans[j];
    		l+=1;
    	}
    		
    	}
		return GeneratorMatrices;
	}
    
    
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FINISH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
   
    
    
	
	 
	 
	
		 
		 
		 
		 
		 
		 
	
		 
		 
		 
		 
		 
		 
		 
		 
		 
	 
	  
	/**
	 * Returns the generator matrices in the natural format, which is 
	 * a 3-dimensional array of shape dim x numRows x numCols.
	 * This format uses an int for each bit, so this is slow and takes much memory!
	 */
	public int[][][] generatorMatricesToStandardFormat(){
		int r, c, j;	// Row r, column c, dimension j.
		int[][][] standardMatrices = new int[dim][numRows][numCols];
		for (j = 0; j < dim; j++) {
			for (c = 0; c < numCols; c++) {
				int column = genMat[j * numCols + c];  // This column as an integer.
				column >>= outDigits - numRows;		   // outDigits (or w) is used here.
				for (r = numRows - 1; r >= 0; r--) {
					standardMatrices[j][r][c] = (column & 1);
					column >>= 1;
				}
			}
		}
		return standardMatrices;
	}

	/**
	 * Sets the generator matrices from matrices in the standard format.
	 * Standard format is a 3-dimensional array of shape dimension x number of rows x number of columns.
	 */
	public void generatorMatricesFromStandardFormat(int[][][] matrices){
		assert matrices.length == dim;
		assert matrices.length > 0;
		assert matrices[0].length == numRows;
		assert matrices[0].length > 0;
		assert matrices[0][0].length == numCols;

		genMat = new int[dim * numCols];
		int r, c, j;	// Row r, column c, dimension j.
		for(j = 0; j < dim; ++j){
			for(r = 0; r < numRows; ++r){
				for(c = 0; c < numCols; ++c){
					if (matrices[j][r][c] > 0)
						genMat[j * numCols + c] += (1 << (outDigits - 1 - r));
				}
			}
		}
	}
	 
	/**
	 * Prints the generator matrices as bit matrices in standard form for dimensions 1 to @f$s@f$.
	 * Each matrix has @f$r@f$ rows and @f$k@f$ columns.
	 */
	public void printGeneratorMatrices(int s) {
		int r, c, j;	// Row r, column c, dimension j.
		int[][][] standardMatrices = generatorMatricesToStandardFormat();
		for (j = 0; j < s; j++) {
			System.out.println("dim = " + (j + 1) + PrintfFormat.NEWLINE);
			for (r = 0; r < numRows; r++) {
				StringBuffer sb = new StringBuffer();
				for (c = 0; c < numCols; c++) {
					sb.append(standardMatrices[j][r][c]);
				}
				System.out.println(sb);
			}
			System.out.println("----------------------------------");
		}
	}

	/**
	 * Prints the generator matrices transposed in the form of integers for dimensions 1 to @f$s@f$.
	 * Each integer corresponds to one column of bits.
	 */
	public void printGeneratorMatricesTrans(int s) {
		// column c, dimension j.
		for (int j = 0; j < s; j++) {
			System.out.println("dim = " + (j + 1) + PrintfFormat.NEWLINE);
			for (int c = 0; c < numCols; c++)
				System.out.println(genMat[j * numCols + c]);
			System.out.println("----------------------------------");
		}
	}

	/**
	 * Returns a copy of the generator matrices transposed in the form of integers for dimensions 1 to @f$s@f$.
	 * Each integer corresponds to one column of bits.
	 */
	public int[] getGeneratorMatricesTrans() {
		return Arrays.copyOf(genMat, genMat.length);
	}

	public double getCoordinate(int i, int j) {
		int res;
		int pos = 0;
		int grayCode = i ^ (i >> 1);

		if (digitalShift == null)
			res = 0;
		else
			res = digitalShift[j];
		while ((grayCode >> pos) != 0) {
			if (((grayCode >> pos) & 1) != 0)
				res ^= genMat[j * numCols + pos];
			pos++;
		}
		if (digitalShift != null)
			return res * normFactor + EpsilonHalf;
		else
			return res * normFactor;
	}

	public double getCoordinateNoGray(int i, int j) {
		int res;
		if (digitalShift == null)
			res = 0;
		else
			res = digitalShift[j];
		int pos = 0;              // Position of the bit that is examined.
		// Add (xor) the columns of C_j for which the corresponding bit of i is 1.
		// Least significant bit of i goes with first column of C_j.
		// The first row of C_j is for the most significant bit of output, u_{i,j,1}.
		while ((i >> pos) != 0) {
			if ((((i >> pos) & 1) != 0) && (pos < numCols))
				res ^= genMat[j * numCols + pos];
			pos++;
		}
		if (digitalShift != null)
			// This EpsilonHalf must be explained in the doc. 
			// It prevents the output of 0.
			// normFactor is 2^{-outDigits}.  
			return res * normFactor + EpsilonHalf;  
		else
			return res * normFactor;
	}

	/**
	 * Returns a `DigitalNetBase2Iterator` which enumerates the points using a Gray code.
	 */
	public PointSetIterator iterator() {
		return new DigitalNetBase2Iterator();
	}

	/**
	 * This iterator does not use the Gray code. Thus the points are enumerated in the order of
	 * their first coordinate before randomization.
	 */
	public PointSetIterator iteratorNoGray() {
		return new DigitalNetBase2IteratorNoGray();
	}

	public String toString() {
		StringBuffer sb = new StringBuffer("DigitalNetBase2: ");
		sb.append(super.toString());
		return sb.toString();
	}

	public void clearRandomShift() {
		super.clearRandomShift();
		digitalShift = null;
	}

	public void addRandomShift(int d1, int d2, RandomStream stream) {
		if (null == stream)
			throw new IllegalArgumentException(
			        PrintfFormat.NEWLINE + "   Calling addRandomShift with null stream");
		if (0 == d2)
			d2 = Math.max(1, dim);
		if (digitalShift == null) {
			digitalShift = new int[d2];
			capacityShift = d2;
		} else if (d2 > capacityShift) {
			int d3 = Math.max(4, capacityShift);
			while (d2 > d3)
				d3 *= 2;
			int[] temp = new int[d3];
			capacityShift = d3;
			for (int i = 0; i < d1; i++)
				temp[i] = digitalShift[i];
			digitalShift = temp;
		}
		int maxj;
		if (outDigits < 31)      // outDigit (= w) is used here for the shift!
			maxj = (1 << outDigits) - 1;
		else
			maxj = 2147483647;
		for (int i = d1; i < d2; i++)
			digitalShift[i] = stream.nextInt(0, maxj);
		dimShift = d2;
		shiftStream = stream;
	}

	public void addRandomShift(RandomStream stream) {
		addRandomShift(0, dim, stream);
	}

	// Left-multiplies lower-triangular matrix Mj by original C_j,
	// where original C_j is in originalMat and result is in genMat.
	// Mj[d] is assumed to contain the d-th subdiagonal of matrix Mj,
	// for d=0,...,w-1. Each subdiagonal is represented as a
	// w-bit integer, whose most significant bits are those on the
	// diagonal. For example, for d=w-3, the subdiagonal has 3 bits,
	// say b1, b2, b3, and is represented by the integer
	// Mj[w-3] = b1 * 2^{w-1} + b2 * 2^{w-2} + b3 * b^{w-3}.
	//
	private void leftMultiplyMat(int j, int[] Mj) {
		int c, d, col;       // Dimension j, column c for new C_j.
		for (c = 0; c < numCols; c++) {
			col = 0;
			for (d = 0; d < outDigits; d++)
				// Multiply subdiagonal d of M_j by column c of C_j, and xor.
				col ^= (Mj[d] & originalMat[j * numCols + c]) >> d;
			genMat[j * numCols + c] = col;
		}
	}

	/*
	 * // Left-multiplies lower-triangular matrix Mj by original C_j, // where original C_j is in
	 * originalMat and result is in genMat. private void leftMultiplyMat (int j, int[] Mj) { int c,
	 * l, i, prod; // Dimension j, column c for new C_j. int numOnes; int col; // Will be column c
	 * of genMat. for (c = 0; c < numCols; c++) { col = 0; for (l = 0; l < outDigits; l++) { //
	 * Multiply row l of M_j by column c of C_j. prod = Mj[l] & originalMat[j* numCols + c]; numOnes
	 * = 0; // Counts the number of ones in prod, modulo 2. for (i = 0; i < outDigits; i++) numOnes
	 * += (1 & prod >> i); // Put a 1 in column c, row l, of C_j if numOnes is odd. col += ((numOnes
	 * & 1) << (outDigits-l-1)); } genMat[j * numCols + c] = col; } }
	 */

	// Right-multiplies upper-triangular matrix Mj by original C_j,
	// where original C_j is in originalMat and result is in genMat.
	// Mj[d] is assumed to contain the d-th column of matrix Mj,
	// for d=0,...,w-1. Each column is represented as a w-bit integer,
	// whose most significant bits are those at index 0.
	// For example, for d=2, the column has 3 bits, (the others are 0
	// since under the diagonal) say b1, b2, b3, and is represented by
	// the integer Mj[2] = b1 * 2^{w-1} + b2 * 2^{w-2} + b3 * b^{w-3}.
	//
	private void rightMultiplyMat(int j, int[] Mj) {
		int c, r, col;       // Dimension j, column c for new C_j.
		int mask;            // Bit of column Mj[c]

		for (c = 0; c < numCols; c++) {
			mask = 1 << outDigits - 1;
			col = originalMat[j * numCols + c];
			for (r = 0; r < c; r++) {
				// If bit (outDigits - 1 - r) of Mj[c] is 1, add column r
				if ((Mj[c] & mask) != 0)
					col ^= originalMat[j * numCols + r];
				mask >>= 1;
			}
			genMat[j * numCols + c] = col;
		}
	}

	/*
	 * // Right-multiplies original C_j by upper-triangular matrix Mj, // where original C_j is in
	 * originalMat and result is in genMat. private void rightMultiplyMat (int j, int[] Mj) { int c,
	 * l, i, mask; // Dimension j, column c for new C_j. int numOnes; int col; // Will be column c
	 * of genMat. boolean bool1, bool2; for (c = 0; c < numCols; c++) { col = 0; for (l = 0; l <
	 * outDigits; l++) { // Multiply row l of C_j by column c of Mj. mask = (1 << outDigits-l-1); //
	 * ??? // xor = originalMat[j* numCols + l] & Mj[c]; numOnes = 0; // Counts the number of ones
	 * in xor, modulo 2. for (i = 0; i < numCols; i++) { bool1 = (mask & originalMat[j * numCols +
	 * i]) != 0; bool2 = ((1 << (outDigits-i-1)) & Mj[i]) != 0; if (bool1 & bool2) numOnes++; } //
	 * Put a 1 in column c, row l, of C_j if numOnes is odd. col += ((numOnes & 1) <<
	 * (outDigits-l-1)); } genMat[j * numCols + c] = col; } }
	 */

	public void leftMatrixScramble(RandomStream stream) {
		int j, d;  // dimension j, subdiagonal d.
		final int allOnes = (1 << outDigits) - 1;    // outDigits ones.

		// If genMat contains the original gen. matrices, copy to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols];
		}
		// Constructs the lower-triangular scrambling matrices M_j, w by w.
		// scrambleMat[j][l] contains row l in a single integer (binary repres.)
		scrambleMat = new int[dim][outDigits];
		for (j = 0; j < dim; j++) {
			scrambleMat[j][0] = allOnes;
			for (d = 1; d < outDigits; d++)
				scrambleMat[j][d] = (stream.nextInt(0, allOnes >> d)) << d;
		}
		
		// Multiply M_j by the generator matrix C_j for each j.
		for (j = 0; j < dim; j++)
			leftMultiplyMat(j, scrambleMat[j]);
	}

	public void iBinomialMatrixScramble(RandomStream stream) {
		int j, d;     // Dimension j, subdiagonal d of M_j.
		final int allOnes = (1 << outDigits) - 1;    // outDigits ones.
		int lastRow;  // Last row of M_j: w-1 random bits followed by 1.

		// If genMat is original generator matrices, copy it to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols];
		}

		// Constructs the lower-triangular scrambling matrices M_j, w by w.
		// scrambleMat[j][l] contains row l of M_j.
		int[][] scrambleMat = new int[dim][outDigits];
		for (j = 0; j < dim; j++) {
			scrambleMat[j][0] = allOnes;
			lastRow = stream.nextInt(0, allOnes) | 1;
			for (d = 1; d < outDigits; d++)
				// Subdiagonal d contains either all ones or all zeros.
				if (((1 << d) & lastRow) == 0)
					scrambleMat[j][d] = 0;
				else
					scrambleMat[j][d] = (allOnes >> d) << d;
		}
		for (j = 0; j < dim; j++)
			leftMultiplyMat(j, scrambleMat[j]);
		// leftMultiplyMat (scrambleMat);
	}

	public void stripedMatrixScramble(RandomStream stream) {
		int j, d;  // dimension j, subdiagonal d of M_j.

		// If genMat is original generator matrices, copy it to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols];
		}
		// Constructs the lower-triangular scrambling matrix M, w by w,
		// filled with 1's. scrambleMat[d] contains subdiagonal d of M.
		int[] scrambleMat = new int[outDigits];
		final int allOnes = (1 << outDigits) - 1;    // outDigits ones.
		for (d = 0; d < outDigits; d++)
			scrambleMat[d] = (allOnes >> d) << d;
		for (j = 0; j < dim; j++)
			leftMultiplyMat(j, scrambleMat);
	}

	/*
	 * public void leftMatrixScramble (RandomStream stream) { int j, l; // dimension j, row l. int
	 * boundInt;
	 * 
	 * // If genMat contains the original gen. matrices, copy to originalMat. if (originalMat ==
	 * null) { originalMat = genMat; genMat = new int[dim * numCols]; }
	 * 
	 * // Constructs the lower-triangular scrambling matrices M_j, w by w. // scrambleMat[j][l]
	 * contains row l in a single integer (binary repres.) int[][] scrambleMat = new
	 * int[dim][outDigits]; for (j = 0 ; j < dim; j++) { boundInt = 0; for (l = 0; l < outDigits;
	 * l++) { boundInt += (1 << l); // Integer repres. by string of l+1 ones. scrambleMat[j][l] =
	 * (stream.nextInt (0, boundInt) | 1) << (outDigits-l-1); } }
	 * 
	 * // Multiply M_j by the generator matrix C_j for each j. for (j = 0; j < dim; j++)
	 * leftMultiplyMat (j, scrambleMat[j]); }
	 * 
	 * public void iBinomialMatrixScramble (RandomStream stream) { int j, l; // dimension j, row l
	 * of M_j. int allOnes;
	 * 
	 * // If genMat is original generator matrices, copy it to originalMat. if (originalMat == null)
	 * { originalMat = genMat; genMat = new int[dim * numCols]; }
	 * 
	 * // Constructs the lower-triangular scrambling matrices M_j, w by w. // scrambleMat[j][l]
	 * contains row l of M_j. int[][] scrambleMat = new int[dim][outDigits]; for (j = 0 ; j < dim;
	 * j++) { allOnes = ~0 >> (32 - outDigits); // outDigits ones. scrambleMat[j][outDigits-1] =
	 * stream.nextInt (0, allOnes) | 1; for (l = outDigits - 2; l >= 0; l--) scrambleMat[j][l] =
	 * scrambleMat[j][l+1] << 1; }
	 * 
	 * for (j = 0; j < dim; j++) leftMultiplyMat (j, scrambleMat[j]); // leftMultiplyMat
	 * (scrambleMat); }
	 * 
	 * public void stripedMatrixScramble (RandomStream stream) { int j, l; // dimension j, row l of
	 * M_j. int allOnes;
	 * 
	 * // If genMat is original generator matrices, copy it to originalMat. if (originalMat == null)
	 * { originalMat = genMat; genMat = new int[dim * numCols]; }
	 * 
	 * // Constructs the lower-triangular scrambling matrix M, w by w, // filled with 1's.
	 * scrambleMat[l] contains row l of M. int[] scrambleMat = new int[outDigits]; allOnes = ~0 >>
	 * (32 - outDigits); // outDigits ones. for (l = 0; l < outDigits; l++) scrambleMat[l] =
	 * (allOnes << (outDigits - 1 - l)) & allOnes; for (j = 0; j < dim; j++) leftMultiplyMat (j,
	 * scrambleMat); }
	 */

	public void rightMatrixScramble(RandomStream stream) {
		int j, c;     // Dimension j, column c for new C_j.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols];
		}
		// Generate an upper triangular matrix for Faure-Tezuka right-scramble.
		// scrambleMat[c] contains column c of M.
		int[] scrambleMat = new int[outDigits];
		int boundInt = 0;
		for (c = 0; c < numCols; c++) {
			boundInt += (1 << c); // Integer repres. by string of c+1 ones.
			scrambleMat[c] = (1 | stream.nextInt(0, boundInt)) << (outDigits - c - 1);
		}
		// Right-multiply the generator matrices by the scrambling matrix.
		for (j = 0; j < dim; j++)
			rightMultiplyMat(j, scrambleMat);
	}

	/**
	 * Generate a vector of `numBits <= 31` random bits using the random stream `stream`.
	 */
	private int randomBitVector(RandomStream stream, int numBits) {
		if (numBits < 1)
			throw new IllegalArgumentException("numBits must be >= 1");
		if (numBits > 31)
			throw new IllegalArgumentException("numBits must be <= 31");
		int maxj;
		if (numBits < 31)
			maxj = (1 << numBits) - 1;
		else
			maxj = 2147483647;
		return stream.nextInt(0, maxj) << (31 - numBits);
	}

	/**
	 * Same as @link nestedUniformScramble(RandomStream,double[][],int)
	 * nestedUniformScramble(stream, output, 0) @endlink.
	 */
	public void nestedUniformScramble(RandomStream stream, double[][] output) {
		nestedUniformScramble(stream, output, 0);
	}

	/**
	 * Applies Owen's nested uniform scrambling to the digital net and returns the scrambled points 
	 * in the two-dimensional array `output`.  The points are computed by using the generating matrices,
	 * and scrambled right away.  Only the first `numBits` bits are scrambled. 
	 * This function does not modify the @ref DigitalNetBase2 object. In particular, it does
	 * not change the generating matrices stored in the object. All points are
	 * randomized at once to avoid storing the permutations.
	 *
	 * The implementation is an adaptation of that found in [SAMPLE
	 * PACKage](http://www.uni-kl.de/AG-Heinrich/SamplePack.html) by Thomas Kollig and Alexander
	 * Keller.
	 *
	 * @param stream
	 *            Random stream used to randomize the bits.
	 * @param output
	 *            Output array that will store the randomized points. The size of its first
	 *            dimension must be getNumPoints() and the size of its second dimension must be
	 *            getDimension().
	 * @param numBits
	 *            Number of output bits to scramble. It can be smaller than, equal to or larger than
	 *            DigitalNet.outDigits. If this parameter is zero, `outDigits` bits will be scrambled.
	 */
	public void nestedUniformScramble(RandomStream stream, double[][] output, int numBits) {
		assert output.length == numPoints;
		assert output.length > 0;
		assert output[0].length == dim;

		int[][] int_output = new int[numPoints][dim];
		nestedUniformScramble(stream, int_output, numBits);
		for (int j = 0; j < dim; ++j) {
			for (int i = 1; i < numPoints; i++) {
				output[i][j] = int_output[i][j] * normFactor + EpsilonHalf;
			}
		}
	}

	/**
	 * Same as @link nestedUniformScramble(RandomStream,double[][],int)@endlink, but it
	 * returns the points as integers between 0 and 2^outDigits - 1, instead of doubles.
	 */
	public void nestedUniformScramble(RandomStream stream, int[][] output, int numBits) {
		assert output.length == numPoints;
		assert output.length > 0;
		assert output[0].length == dim;

		if (numBits == 0)
			numBits = outDigits;

		int[] poslist = new int[2 * numPoints];
		int[] bvlist = new int[2 * numPoints];
		int[] counts = new int[256];
		int[] binpos = new int[256];

		for (int j = 0; j < dim; ++j) {
			bvlist[0] = 0;
			poslist[0] = 0;
			for (int i = 1; i < numPoints; i++) {
				// We use Gray code order (could be optional).
				// We could have used a point set iterator here, but the iterator computes all
				// coordinates at once and we need only one at a time.
				int pos = 0;
				int bv = 1;
				while ((i & bv) == 0) {
					pos++;
					bv <<= 1;
				}
				bvlist[i] = bvlist[i - 1] ^ genMat[j * numCols + pos];
				poslist[i] = i;
			}
			for (int b = 0; b < 4; b++) {
				for (int i = 0; i < 256; i++)
					counts[i] = 0;
				int m = (b % 2) * numPoints;
				int bb = 8 * b;
				int bv = 0xff << bb;
				for (int i = 0; i < numPoints; i++)
					counts[(bvlist[m + i] & bv) >>> bb]++;
				binpos[0] = (1 - b % 2) * numPoints;
				for (int i = 0; i < 255; i++)
					binpos[i + 1] = binpos[i] + counts[i];
				for (int i = 0; i < numPoints; i++) {
					int pos = (bvlist[m + i] & bv) >>> bb;
					int k = binpos[pos]++;
					bvlist[k] = bvlist[m + i];
					poslist[k] = poslist[m + i];
				}
			}
			int bv = randomBitVector(stream, numBits);
			output[poslist[0]][j] = bvlist[0] ^ bv;
			for (int i = 1; i < numPoints; i++) {
				int bv2 = bvlist[i - 1];
				bv2 ^= bvlist[i];
				bv2 = randomBitVector(stream, numBits)
				        & ((int) (1 << (int) Num.log2((double) bv2)) - 1);
				bv ^= bv2;
				output[poslist[i]][j] = bvlist[i] ^ bv;
			}

		}
	}

	// -----------------------------------------------------------------------
	private void ScrambleError(String method) {
		throw new UnsupportedOperationException(
		        PrintfFormat.NEWLINE + "  " + method + " is not implemented for a DigitalNetBase2");
	}

	public void leftMatrixScrambleDiag(RandomStream stream) {
		ScrambleError("leftMatrixScrambleDiag");
	}

	public void leftMatrixScrambleFaurePermut(RandomStream stream, int sb) {
		ScrambleError("leftMatrixScrambleFaurePermut");
	}

	public void leftMatrixScrambleFaurePermutDiag(RandomStream stream, int sb) {
		ScrambleError("leftMatrixScrambleFaurePermutDiag");
	}

	public void leftMatrixScrambleFaurePermutAll(RandomStream stream, int sb) {
		ScrambleError("leftMatrixScrambleFaurePermutAll");
	}

	public void iBinomialMatrixScrambleFaurePermut(RandomStream stream, int sb) {
		ScrambleError("iBinomialMatrixScrambleFaurePermut");
	}

	public void iBinomialMatrixScrambleFaurePermutDiag(RandomStream stream, int sb) {
		ScrambleError("iBinomialMatrixScrambleFaurePermutDiag");
	}

	public void iBinomialMatrixScrambleFaurePermutAll(RandomStream stream, int sb) {
		ScrambleError("iBinomialMatrixScrambleFaurePermutAll");
	}

	public void stripedMatrixScrambleFaurePermutAll(RandomStream stream, int sb) {
		ScrambleError("stripedMatrixScrambleFaurePermutAll");
	}


	/**
	 * Interlaces the points from a digital net.
	 *
	 * This is useful for interlacing after NUS, since NUS returns the points. For other use
	 * cases, @link matrixInterlace()@endlink should be preferred.
	 *
	 * @param points
	 *            Array that stores the non-interlaced points as integers. The size of its first
	 *            dimension must be getNumPoints() and the size of its second dimension must be
	 *            getDimension().
	 * @param interlacedPoints
	 *            Output array that will store the interlaced points. The size of its first
	 *            dimension must be getNumPoints() and the size of its second dimension must be
	 *            getDimension() / getInterlacing().
	 */
	public void outputInterlace(int[][] points, double[][] interlacedPoints) {
		assert points.length == numPoints;
		assert points.length > 0;
		assert points[0].length == dim;

		assert interlacedPoints.length == numPoints;
		assert interlacedPoints.length > 0;
		assert interlacedPoints[0].length * interlacing == dim;

		double longNormFactor = 1 / (Math.pow(2, 63));

		for (int i = 0; i < numPoints; i++) {
			for (int j = 0; j < dim / interlacing; ++j){
				long result = 0;
				for (int idx = 0; idx < interlacing; ++idx){
					int coord = j * interlacing + idx;
					int interlaced_pos = 62 - idx;
					int original_pos = 30;
					int mask = 1 << original_pos;
					while (interlaced_pos >= 0 && original_pos >= 0){
						assert interlaced_pos >= original_pos;
						result += (((long) (points[i][coord] & mask)) << (interlaced_pos - original_pos));
						interlaced_pos -= interlacing;
						original_pos -= 1;
						mask >>= 1;
					}
				}
				interlacedPoints[i][j] = result * longNormFactor + EpsilonHalf;
			}
		}
	}

	/**
	 * Interlaces the matrices from a digital net.
	 *
	 * This function returns a new digital net, whose dimension equals getDimension() / getInterlacing(),
	 * and whose generating matrices are interlaced.
	 */
	public DigitalNetBase2 matrixInterlace() {
		DigitalNetBase2 result = new DigitalNetBase2();
		result.dim = dim / interlacing;
		result.interlacing = 1;
		result.numCols = numCols;
		result.numRows = Math.min(MAXBITS, numRows*interlacing);
		result.numPoints = numPoints;
		result.outDigits = outDigits;
		result.normFactor = normFactor;

		int[][][] nonInterlacedMatrices = generatorMatricesToStandardFormat();
		int[][][] interlacedMatrices = new int[result.dim][result.numRows][result.numCols];
		int r, j;	// Row r, dimension j.
		for(j = 0; j < result.dim; ++j){
			for(r = 0; r < result.numRows; ++r){
				interlacedMatrices[j][r] = nonInterlacedMatrices[j * interlacing + r % interlacing][r / interlacing];
			}
		}
		result.generatorMatricesFromStandardFormat(interlacedMatrices);
		
		return result;
		
	}

	
	// *******************************************************************
	protected class DigitalNetBase2Iterator extends DigitalNetIterator {

		// Coordinates of the current point stored (cached) as integers.
		// Initially contains zeros, because first point is always zero.
		// Incorporates the random shift, except for the first point.
		// There is one more dimension for the points because of the
		// shift iterators children of DigitalNetBase2Iterator.
		// dimS = dim, except for the shift iterator children where
		// dimS = dim + 1.
		protected int dimS;

		public DigitalNetBase2Iterator() {
			super();
			EpsilonHalf = 0.5 / Num.TWOEXP[outDigits];
			cachedCurPoint = new int[dim + 1];
			dimS = dim;
			init2();
		}

		public void init() {   // This method is necessary to overload
		}                      // the init() of DigitalNetIterator

		public void init2() { // See constructor
			resetCurPointIndex();
		}

		// We want to avoid generating 0 or 1
		public double nextDouble() {
			return nextCoordinate();
		}

		public double nextCoordinate() {
			if (curPointIndex >= numPoints || curCoordIndex >= dimS)
				outOfBounds();
			if (digitalShift == null)
				return cachedCurPoint[curCoordIndex++] * normFactor;
			else
				return cachedCurPoint[curCoordIndex++] * normFactor + EpsilonHalf;
		}

		protected void addShiftToCache() {
			if (digitalShift == null)
				for (int j = 0; j < dim; j++)
					cachedCurPoint[j] = 0;
			else {
				if (dimShift < dimS)
					addRandomShift(dimShift, dimS, shiftStream);
				for (int j = 0; j < dim; j++)
					cachedCurPoint[j] = digitalShift[j];
			}
		}

		public void resetCurPointIndex() {
			addShiftToCache();
			curPointIndex = 0;
			curCoordIndex = 0;
		}

		public void setCurPointIndex(int i) {
			if (i == 0) {
				resetCurPointIndex();
				return;
			}
			// Out of order computation, must recompute the cached current
			// point from scratch.
			curPointIndex = i;
			curCoordIndex = 0;
			addShiftToCache();

			int j;
			int grayCode = i ^ (i >> 1);
			int pos = 0;      // Position of the bit that is examined.
			while ((grayCode >> pos) != 0) {
				if (((grayCode >> pos) & 1) != 0)
					for (j = 0; j < dim; j++)
						cachedCurPoint[j] ^= genMat[j * numCols + pos];
				pos++;
			}
		}

		public int resetToNextPoint() {
			int pos = 0;  // Will be position of change in Gray code,
			              // = pos. of first 0 in binary code of point index.
			while (((curPointIndex >> pos) & 1) != 0)
				pos++;
			if (pos < numCols) {
				for (int j = 0; j < dim; j++)
					cachedCurPoint[j] ^= genMat[j * numCols + pos];
			}
			curCoordIndex = 0;
			return ++curPointIndex;
		}

		public int nextPoint(double p[], int d) {
			if (curPointIndex >= numPoints || d > dimS)
				outOfBounds();
			if (digitalShift == null) {
				for (int j = 0; j < d; j++)
					p[j] = cachedCurPoint[j] * normFactor;
			} else {
				for (int j = 0; j < d; j++)
					p[j] = cachedCurPoint[j] * normFactor + EpsilonHalf;
			}
			return resetToNextPoint();
		}
	}

	
	// *******************************************************************
	protected class DigitalNetBase2IteratorNoGray extends DigitalNetBase2Iterator {

		// Same as DigitalNetBase2Iterator,
		// except that the Gray code is not used.

		public DigitalNetBase2IteratorNoGray() {
			super();
		}

		public void setCurPointIndex(int i) {
			if (i == 0) {
				resetCurPointIndex();
				return;
			}
			// Out of order computation, must recompute the cached current
			// point from scratch.
			curPointIndex = i;
			curCoordIndex = 0;
			addShiftToCache();
			int pos = 0;      // Position of the bit that is examined.
			while ((i >> pos) != 0) {
				if ((((i >> pos) & 1) != 0) && (pos < numCols)) {
					for (int j = 0; j < dim; j++)
						cachedCurPoint[j] ^= genMat[j * numCols + pos];
				}
				pos++;
			}
		}

		public int resetToNextPoint() {
			// Contains the bits of i that changed.
			if (curPointIndex + 1 >= numPoints)
				return ++curPointIndex;
			int diff = curPointIndex ^ (curPointIndex + 1);
			int pos = 0;      // Position of the bit that is examined.
			while ((diff >> pos) != 0) {
				if ((((diff >> pos) & 1) != 0) && (pos < numCols)) {
					for (int j = 0; j < dim; j++)
						cachedCurPoint[j] ^= genMat[j * numCols + pos];
				}
				pos++;
			}
			curCoordIndex = 0;
			return ++curPointIndex;
		}

	}
}
