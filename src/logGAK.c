// Adapted from the code available at http://marcocuturi.net/GA.html

/*
 ***** BEGIN LICENSE BLOCK *****
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * The Original Code is Global Alignment Kernel, (C) 2010, Marco Cuturi
 *
 * The Initial Developers of the Original Code is
 *
 * Marco Cuturi   mcuturi@i.kyoto-u.ac.jp
 *
 * Portions created by the Initial Developers are
 * Copyright (C) 2011 the Initial Developers. All Rights Reserved.
 *
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 ***** END LICENSE BLOCK *****
 *
 * REVISIONS:
 * This version is functionally equivalent to v1.03, adapted to be called from R. The 'logs'
 * variable is allocated in R (see its purpose below). The 'LOGP' macro was changed to an in-line
 * function.
 *
 * Previous versions:
 * v1.03 Added log1p function for windows platforms, September 12th 2011.
 * v1.02 Changed some C syntax that was not compiled properly on Windows platforms, June 8th
 * v1.01 of Global Alignment Kernel, May 12th 2011 (updated comments fields)
 * v1.0 of Global Alignment Kernel, March 25th 2011.
 */

#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "dtwclust.h"

// Useful constants
#define LOG0 -10000          // log(0)

// LOGP
double inline LOGP(double const x, double const y) __attribute__((always_inline));

double inline LOGP(double const x, double const y) {
    return (x > y) ? x + log1p(exp(y - x)) : y + log1p(exp(x - y));
}

// global alignment kernel
double logGAK_c(double const *seq1 , double const *seq2,
                int const nX, int const nY, int const dim,
                double const sigma, int const triangular,
                double *logs)
{
    /*
     * Implementation of the (Triangular) global alignment kernel.
     *
     * See details about the wrapper function below for more information on the inputs that need
     * to be called from R.
     *
     * seq1 is a first sequence represented as a matrix of real elements.
     *   Each line i corresponds to the vector of observations at time i.
     * seq2 is the second sequence formatted in the same way.
     * nX, nY and dim provide the number of lines of seq1 and seq2.
     * sigma stands for the bandwidth of the \phi_\sigma distance used kernel
     * triangular is a parameter which parameterizes the triangular kernel
     */

    int i, j, ii, cur, old, frompos1, frompos2, frompos3;
    double gram, Sig, aux, sum;

    int curpos = 0;                 // R compilation in windows gives warning about uninitialization...
    int trimax = (nX > nY) ? nX - 1 : nY - 1; // Maximum of abs(i-j) when 1<=i<=nX and 1<=j<=nY
    int cl = trimax + 2;            // length of a column for the dynamic programming

    // logs is the array that will stores two successive columns of the (nX+1) x (nY+1)
    // table used to compute the final kernel value, as well as the triangular coefficients
    // in the third 'column'

    // initialize
    ii = (trimax < triangular) ? trimax + 1 : triangular;
    aux = triangular ? LOG0 : 0;

    for (i = 0; i <= trimax; i++) {
        if (triangular > 0 && i < ii)
            logs[i + 2*cl] = log(1 - i / triangular);
        else
            logs[i + 2*cl] = aux;
    }

    Sig = -1 / (2 * sigma * sigma);

    /****************************************************/
    /* First iteration : initialization of columns to 0 */
    /****************************************************/
    // The left most column is all zeros...
    for (j = 1; j < cl; j++) {
        logs[j] = LOG0;
    }

    // ... except for the lower-left cell which is initialized with a value of 1, i.e. a log value of 0.
    logs[0] = 0;

    // Cur and Old keep track of which column is the current one and which one is the already computed one.
    cur = 1;      // Indexes [0..cl-1] are used to process the next column
    old = 0;      // Indexes [cl..2*cl-1] were used for column 0

    /************************************************/
    /* Next iterations : processing columns 1 .. nX */
    /************************************************/
    // Main loop to vary the position for i=1..nX
    for (i = 1; i <= nX; i++) {
        // Special update for positions (i=1..nX,j=0)
        curpos = cur * cl;                  // index of the state (i,0)
        logs[curpos] = LOG0;

        // Secondary loop to vary the position for j=1..nY
        for (j = 1; j <= nY; j++) {
            curpos = cur * cl + j;          // index of the state (i,j)

            if (logs[abs(i-j) + 2*cl] > LOG0) {
                frompos1 = old*cl + j;      // index of the state (i-1,j)
                frompos2 = cur*cl + j-1;    // index of the state (i,j-1)
                frompos3 = old*cl + j-1;    // index of the state (i-1,j-1)

                // We first compute the kernel value
                sum = 0;
                for (ii = 0; ii < dim; ii++) {
                    sum += pow(seq1[i-1 + ii*nX] - seq2[j-1 + ii*nY], 2);
                }

                gram = logs[abs(i-j) + 2*cl] + sum*Sig;
                gram -= log(2 - exp(gram));

                // Doing the updates now, in two steps
                aux = LOGP(logs[frompos1], logs[frompos2]);
                logs[curpos] = LOGP(aux, logs[frompos3]) + gram;

            } else {
                logs[curpos] = LOG0;
            }
        }

        // Update the culumn order
        cur = 1 - cur;
        old = 1 - old;
    }

    aux = logs[curpos];

    // Return the logarithm of the Global Alignment Kernel
    return aux;
}

// The gateway function
SEXP logGAK(SEXP x, SEXP y, SEXP nx, SEXP ny, SEXP dim, SEXP sigma, SEXP window, SEXP logs) {
    /*
     * Inputs are, in this order
     * A N1 x d matrix (d-variate time series with N1 observations)
     * A N2 x d matrix (d-variate time series with N2 observations)
     * The length of series N1 and N2, and the number of variables d
     * A sigma > 0 parameter for the Gaussian kernel's width
     * A triangular integer which parameterizes the band.
     *   - when triangular = 0, the triangular kernel is not used, the evaluation is on the sum of all
     *   paths, and the complexity if of the order of N1 x N2
     *   - when triangular > 0, the triangular kernel is used and downweights some of the paths that
     *   lay far from the diagonal.
     * A (max(N1,N2) + 1) x 3 matrix of doubles
     */

    int triangular = asInteger(window);
    int nX = asInteger(nx);
    int nY = asInteger(ny);

    double d;

    // If triangular is smaller than the difference in length of the time series, the kernel is equal to zero,
    // i.e. its log is set to -Inf
    if (triangular > 0 && abs(nX - nY) > triangular)
        d = R_NegInf;
    else
        d = logGAK_c(REAL(x), REAL(y), nX, nY, asInteger(dim), asReal(sigma), triangular, REAL(logs));

    return ScalarReal(d);
}
