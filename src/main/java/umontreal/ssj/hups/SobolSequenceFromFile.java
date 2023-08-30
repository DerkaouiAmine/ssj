/*
 * Class:        SobolSequenceFromFile
 * Description:  read the parameters defining a Sobol sequence
 *               from a file or from a URL address
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

import java.io.*;
import java.net.MalformedURLException;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class permits one to read the parameters that define a Sobol sequence
 * either from a file or from a URL address.
 * The parameters used in building the sequence are those defined in the class
 * @ref umontreal.ssj.hups.SobolSequence.
 * The format of the data files must be as follows:

 * <center>
 *  <tt><div class="SSJ-fbox">
 * <table class="SSJ-table">
 * <tr>
 *   <td colspan="2" class="l">// Any number of comment lines starting with
 * //</td>
 * </tr><tr>
 *   <td class="l">@f$d@f$</td>
 *   <td class="l">// Number of dimensions</td>
 * </tr><tr>
 *   <td class="l">@f$n@f$</td>
 *   <td class="l">// Number of points = @f$2^d@f$</td>
 * </tr><tr>
 *   <td colspan="2" class="l">// dim = 1</td>
 * </tr><tr>
 *   <td class="l">@f$a@f$</td>
 *   <td class="l">// @f$= 2^{30}B_1 + 2^{29}B_2 + \cdots+ 2^{31 - d}B_d@f$</td>
 * </tr><tr>
 *   <td colspan="2" class="l">// dim = 2</td>
 * </tr><tr>
 *   <td class="l">@f$\vdots@f$</td>
 *   <td class="l"></td>
 * </tr><tr>
 *   <td colspan="2" class="l">// dim = @f$d@f$</td>
 * </tr><tr>
 *   <td class="l">@f$a@f$</td>
 *   <td class="l"></td>
 * </tr>
 * </table>
 *  </div> </tt>
 * </center>
 *
 * For each dimension @f$j@f$, there must be a 32-bit integer @f$a@f$ corresponding
 * to the coefficient @f$B_j@f$ of the primitive polynomial @f$B_j(x)@f$ used to
 * compute the direction numbers of the sequence.
 * The values of @f$a@f$ must be input in hexadecimal format starting with
 * `0x' or `0X'. For example,
 * <center>
 * <tt><div class="SSJ-fbox">
 * <table class="SSJ-table">
 * <tr>
 *   <td class="l">// Any number of comment lines</td>
 * </tr><tr>
 *   <td class="l">3</td>
 *   <td class="l">// Number of dimensions = 3</td>
 * </tr><tr>
 *   <td class="l">8</td>
 *   <td class="l">// Number of points = @f$2^3 = 8@f$</td>
 * </tr><tr>
 *   <td class="l">0x0000000D</td>
 *   <td class="l">// @f$2^{30}B_1 + 2^{29}B_2 + 2^{28}B_3@f$</td>
 * </tr><tr>
 *   <td class="l">0x0000001B</td>
 *   <td class="l">// @f$2^{30}B_1 + 2^{29}B_2 + 2^{31}B_3@f$</td>
 * </tr><tr>
 *   <td class="l">0x0000002D</td>
 *   <td class="l">// @f$2^{30}B_1 + 2^{28}B_2 + 2^{31}B_3@f$</td>
 * </tr>
 * </table>
 * </div></tt>
 * </center>
 *

 *
 * The method `toString' returns a string with all the parameters read, as
 * well as the name of the file.
 *
 * The method `toStringDetailed' returns a more detailed string including
 * the generating matrices of the sequence.
 *
 * The method `listDir' lists all the files (or directories) in a specified
 * directory. This may be useful to get the list of parameter files available
 * for digital nets.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class SobolSequenceFromFile extends SobolSequence {
	   private String filename;

	   // s is the effective dimension if > 0, otherwise it is dim
	   private void readData (Reader re, int r1, int s1)
	       throws IOException, NumberFormatException
	   {
	      try {
	         StreamTokenizer st = new StreamTokenizer (re);
	         if (st == null) return;
	         st.eolIsSignificant (false);
	         st.slashSlashComments (true);

	         int i = st.nextToken ();
	         if (i != StreamTokenizer.TT_NUMBER)
	            throw new NumberFormatException();
	         b = (int) st.nval;
	         st.nextToken ();   numCols = (int) st.nval;
	         st.nextToken ();   numRows = (int) st.nval;
	         st.nextToken ();   numPoints = (int) st.nval;
	         st.nextToken ();   dim = (int) st.nval;
	         if (dim < 1) {
	            System.err.println (PrintfFormat.NEWLINE +
	                "SobolSequenceFromFile:   dimension dim <= 0");
	            throw new IllegalArgumentException ("dimension dim <= 0");
	         }
	         if (r1 > numRows)
	            throw new IllegalArgumentException (
	            "SobolSequenceFromFile:   One must have   r1 <= Max num rows");
	         if (s1 > dim) {
	            throw new IllegalArgumentException ("s1 is too large");
	         }
	         if (s1 > 0)
	            dim = s1;
	         if (r1 > 0)
	            numRows = r1;

	         if (b != 2) {
	            System.err.println (
	              "***** SobolSequenceFromFile:    only base 2 allowed");
	            throw new IllegalArgumentException ("only base 2 allowed");
	         }
	         genMat = new int[dim * numCols];
	         for (i = 0; i < dim; i++)
	            for (int c = 0; c < numCols; c++) {
	                st.nextToken ();
	                genMat[i*numCols + c] = (int) st.nval;
	            }

	      } catch (NumberFormatException e) {
	         System.err.println (
	            "   SobolSequenceFromFile:   not a number  " + e);
	         throw e;
	      }
	   }


	    private void maskRows (int r, int w) {
	       // Keep only the r most significant bits and set the others to 0.
	       int mask = (int) ((1L << r) - 1);
	       mask <<= MAXBITS - r;
	       for (int i = 0; i < dim; i++)
	          for (int c = 0; c < numCols; c++) {
	              genMat[i*numCols + c] &= mask;
	              genMat[i*numCols + c] >>= MAXBITS - w;
	          }
	    }

	   /**
	    * Constructs a Sobol Sequence after reading its parameters from
	    * file `filename`. See the documentation in
	    * @ref umontreal.ssj.hups.DigitalNetFromFile. Parameter `w` gives the
	    * number of bits of resolution, `r1` is the number of rows, and `s1`
	    * is the dimension. Restrictions: `s1` must be less than the maximal
	    * dimension, and `r1` less than the maximal number of rows in the data
	    * file. Also `w` @f$\ge@f$ `r1`.
	    *  @param filename     Name of the file to be read
	    *  @param r1           Number of rows for the generating matrices
	    *  @param w            Number of bits of resolution
	    *  @param s1           Number of dimensions
	    */
	   public SobolSequenceFromFile (String filename, int r1, int w, int s1)
	         throws IOException, MalformedURLException
	   {
	      super (r1, w, s1);
	      if (w < r1 || w > MAXBITS)
	         throw new IllegalArgumentException (" Must have numRows <= w <= 31");

	      BufferedReader input;
	      if (filename.startsWith("http:") || filename.startsWith("ftp:"))
	         input = DigitalNetFromFile.openURL(filename);
	      else
	         input = DigitalNetFromFile.openFile(filename);

	      try {
	         readData (input, r1, s1);
	      } catch (NumberFormatException e) {
	         System.err.println (
	            "   SobolSequenceFromFile:   cannot read from   " + filename);
	         throw e;

	      }  catch (IOException e) {
	         System.err.println (
	            "   SobolSequenceFromFile:  cannot read from  " + filename);
	         throw e;
	      }
	      input.close();
	      maskRows (numRows, w);
	      outDigits = w;
	      if (numCols >= MAXBITS)
	         throw new IllegalArgumentException (" Must have numCols < 31");

	      this.filename = filename;
	      int x = (1 << numCols);
	      if (x != numPoints) {
	         System.out.println ("numPoints != 2^k");
	         throw new IllegalArgumentException ("numPoints != 2^k");
	      }
	      // Compute the normalization factors.
	      normFactor = 1.0 / ((double) (1L << (outDigits)));

	  }

	   /**
	    * Same as  {@link #SobolSequenceFromFile()
	    * SobolSequenceFromFile(filename, r, 31, s1)} where `s1` is the
	    * dimension and `r` is given in data file `filename`.
	    *  @param filename     Name of the file to be read
	    *  @param s1           Number of dimensions
	    */
	   public SobolSequenceFromFile(String filename, int s1)
	        throws IOException, MalformedURLException
	   {
	       this (filename, -1, 31, s1);
	   }


	   public String toString() {
	      StringBuffer sb = new StringBuffer ("File:  " + filename  +
	         PrintfFormat.NEWLINE);
	      sb.append (super.toString());
	      return sb.toString();
	   }

	/**
	 * Writes the parameters and the generating matrices of this digital net to a
	 * string. This is useful to check that the file parameters have been read
	 * correctly.
	 */
	public String toStringDetailed() {
	      StringBuffer sb = new StringBuffer (toString() + PrintfFormat.NEWLINE);
	      sb.append ("dim = " + dim + PrintfFormat.NEWLINE);
	      for (int i = 0; i < dim; i++) {
	         sb.append (PrintfFormat.NEWLINE + "// dim = " + (1 + i) +
	              PrintfFormat.NEWLINE);
	         for (int c = 0; c < numCols; c++)
	            sb.append  (genMat[i*numCols + c]  + PrintfFormat.NEWLINE);
	      }
	      sb.append ("--------------------------------" + PrintfFormat.NEWLINE);
	      return sb.toString ();
	   }

	   /**
	    * Lists all files (or directories) in directory `dirname`. Only
	    * relative pathnames should be used. The files are parameter files
	    * used in defining digital nets. For example, calling `listDir("")`
	    * will give the list of the main data directory in SSJ, while calling
	    * `listDir("Edel/OOA2")` will give the list of all files in directory
	    * `Edel/OOA2`.
	    */
	   public static String listDir (String dirname) throws IOException {
	      return DigitalNetFromFile.listDir(dirname);
	   }

}
