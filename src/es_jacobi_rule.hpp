// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                   Gauss quadrature formulas                                               *****
// *****   Original C++ code by J. Burkadrt (people.sc.fsu.edu/~jburkardt/cpp_src/jacobi_rule/jacobi_rule.html)    *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#ifndef ES_JACOBI_RULE_HPP
    #define ES_JACOBI_RULE_HPP
#endif

//****************************************************************************80
// Gamma-function - for particular values
//****************************************************************************80
template<typename fpv>
fpv tgamma(fpv arg) {
    fpv tinyloc = get_eps<fpv>() * 20.0;
    fpv sqrt_pi = sqrt(GetPiNumber<fpv>());
    if(fabs(arg - 0.5) < tinyloc) return sqrt_pi;
    if(fabs(arg - 1.0) < tinyloc) return 1.;
    if(fabs(arg - 1.5) < tinyloc) return 0.5*sqrt_pi;
    if(fabs(arg - 2.0) < tinyloc) return 1.;
    if(fabs(arg - 2.5) < tinyloc) return 0.75*sqrt_pi;
    if(fabs(arg - 3.0) < tinyloc) return 2;
    if(fabs(arg - 3.5) < tinyloc) return 1.875*sqrt_pi;
    crash("tgamma: not implemented");
}
//****************************************************************************80


//****************************************************************************80
// PARCHK checks parameters ALPHA and BETA for classical weight functions.
//****************************************************************************80
template<typename fpv>
void parchk ( int kind, int m, fpv alpha, fpv beta ) {
    ASSERT(kind > 0, "KIND <= 0");
    // Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
    if(3 <= kind) ASSERT(alpha > -1.0, "3 <= KIND and ALPHA <= -1.");
    //  Check BETA for Jacobi.
    if ( kind == 4) ASSERT(beta > -1.0, "KIND == 4 and BETA <= -1.0");
    //  Check ALPHA and BETA for rational.
    if ( kind == 8 ) { fpv tmp = alpha + beta + m + 1.0; ASSERT(0.0 > tmp && tmp > beta, "KIND == 8 but condition on ALPHA and BETA fails."); }
}
//****************************************************************************80


//****************************************************************************80
//
//  Purpose:
//
//    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
//
//  Discussion:
//
//    This routine computes the diagonal AJ and sub-diagonal BJ
//    elements of the order M tridiagonal symmetric Jacobi matrix
//    associated with the polynomials orthogonal with respect to
//    the weight function specified by KIND.
//
//    For weight functions 1-7, M elements are defined in BJ even
//    though only M-1 are needed.  For weight function 8, BJ(M) is
//    set to zero.
//
//    The zero-th moment of the weight function is returned in ZEMU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, int M, the order of the Jacobi matrix.
//    Input, double ALPHA, the value of Alpha, if needed.
//    Input, double BETA, the value of Beta, if needed.
//    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
//        of the Jacobi matrix.
//    Output, double CLASS_MATRIX, the zero-th moment.
//
//****************************************************************************80

template<typename fpv>
fpv class_matrix ( int kind, int m, fpv alpha, fpv beta, fpv* aj, fpv* bj) {
    parchk ( kind, 2 * m - 1, alpha, beta );
    fpv zemu; // int p(x) dx

    switch(kind) {
    case 1: {
        fpv ab = 0.0;
        zemu = 2.0 / ( ab + 1.0 );
        for (int i = 0; i < m; i++ ) aj[i] = 0.0;

        for (int i = 1; i <= m; i++ ) {
            fpv abi = i + ab * ( i % 2 );
            fpv abj = 2 * i + ab;
            bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
        }
        break;
    }
    case 2: {
        zemu = GetPiNumber<fpv>();
        for (int i = 0; i < m; i++ ) aj[i] = 0.0;
        bj[0] =  sqrt ( 0.5 );
        for (int i = 1; i < m; i++ ) bj[i] = 0.5;
        break;
    }
    case 3: {
        fpv ab = alpha * 2.0;
        zemu = pow ( 2.0, ab + 1.0 ) * SQR ( tgamma ( alpha + 1.0 ) ) / tgamma ( ab + 2.0 );
        for (int i = 0; i < m; i++ ) aj[i] = 0.0;
        bj[0] = sqrt ( 1.0 / ( 2.0 * alpha + 3.0 ) );
        for (int i = 2; i <= m; i++ )
            bj[i-1] = sqrt ( i * ( i + ab ) / ( 4.0 * SQR ( i + alpha ) - 1.0 ) );
        break;
    }
    case 4: {
        fpv ab = alpha + beta;
        fpv abi = 2.0 + ab;
        zemu = pow ( 2.0, ab + 1.0 ) * tgamma ( alpha + 1.0 ) * tgamma ( beta + 1.0 ) / tgamma ( abi );
        aj[0] = ( beta - alpha ) / abi;
        bj[0] = sqrt ( 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) / ( ( abi + 1.0 ) * abi * abi ) );
        fpv a2b2 = beta * beta - alpha * alpha;

        for (int i = 2; i <= m; i++ ) {
            abi = 2.0 * i + ab;
            aj[i-1] = a2b2 / ( ( abi - 2.0 ) * abi );
            abi = abi * abi;
            bj[i-1] = sqrt ( 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) / ( ( abi - 1.0 ) * abi ) );
        }
        break;
    }
    case 5: {
        zemu = tgamma ( alpha + 1.0 );
        for (int i = 1; i <= m; i++ ) {
            aj[i-1] = 2.0 * i - 1.0 + alpha;
            bj[i-1] = sqrt ( i * ( i + alpha ) );
        }
        break;
    }
    case 6: {
        zemu = tgamma ( ( alpha + 1.0 ) / 2.0 );
        for (int i = 0; i < m; i++ ) aj[i] = 0.0;
        for (int i = 1; i <= m; i++ )
            bj[i-1] = sqrt ( ( i + alpha * ( i % 2 ) ) / 2.0 );
        break;
    }
    case 7: {
        fpv ab = alpha;
        zemu = 2.0 / ( ab + 1.0 );
        for (int i = 0; i < m; i++ ) aj[i] = 0.0;
        for (int i = 1; i <= m; i++ ) {
            fpv abi = i + ab * ( i % 2 );
            fpv abj = 2 * i + ab;
            bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
        }
        break;
    }
    case 8: {
        fpv ab = alpha + beta;
        zemu = tgamma ( alpha + 1.0 ) * tgamma ( - ( ab + 1.0 ) ) / tgamma ( - beta );
        fpv apone = alpha + 1.0;
        fpv aba = ab * apone;
        aj[0] = - apone / ( ab + 2.0 );
        bj[0] = - aj[0] * ( beta + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
        for (int i = 2; i <= m; i++ ) {
            fpv abti = ab + 2.0 * i;
            aj[i-1] = aba + 2.0 * ( ab + i ) * ( i - 1 );
            aj[i-1] = - aj[i-1] / abti / ( abti - 2.0 );
        }

        for (int i = 2; i <= m - 1; i++ ) {
            fpv abti = ab + 2.0 * i;
            bj[i-1] = i * ( alpha + i ) / ( abti - 1.0 ) * ( beta + i ) / ( abti * abti ) * ( ab + i ) / ( abti + 1.0 );
        }
        bj[m-1] = 0.0;
        for (int i = 0; i < m; i++ ) bj[i] =  sqrt ( bj[i] );
        break;
    }
    default: crash("Class-matrix: wrong kind of quadrature");
    }
    return zemu;
}
//****************************************************************************80


//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine.
//
//    It has been modified to produce the product Q' * Z, where Z is an input
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.
//    The changes consist (essentially) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
//****************************************************************************80
template<typename fpv>
void imtqlx ( int n, fpv* d, fpv* e, fpv* z ) {
    const fpv prec = get_eps<fpv>();
    const int itn = 30;

    if ( n == 1 ) return;

    e[n-1] = 0.0;

    int m = -1;
    for (int l = 1; l <= n; l++ ) {
        int j = 0;
        for ( ; ; ) {
            for (m = l; m <= n; m++ ) {
                if ( m == n ) break;
                if ( fabs ( e[m-1] ) <= prec * ( fabs ( d[m-1] ) + fabs ( d[m] ) ) ) break;
            }
            if ( m == l ) break;

            if ( itn <= j ) crash("IMTQLX - Fatal error! Iteration limit exceeded");
            j++;
            fpv g = ( d[l] - d[l-1] ) / ( 2.0 * e[l-1] );
            fpv r =  sqrt ( g * g + 1.0 );
            g = d[m-1] - d[l-1] + e[l-1] / ( g + fabs ( r ) * (g<0.0 ? -1.0 : 1.0) );
            fpv s = 1.0;
            fpv c = 1.0;
            fpv p = 0.0;
            int mml = m - l;

            for (int ii = 1; ii <= mml; ii++ ) {
                int i = m - ii;
                fpv f = s * e[i-1];
                fpv b = c * e[i-1];

                if ( fabs ( g ) <= fabs ( f ) ) {
                    c = g / f;
                    r =  sqrt ( c * c + 1.0 );
                    e[i] = f * r;
                    s = 1.0 / r;
                    c = c * s;
                }
                else {
                    s = f / g;
                    r =  sqrt ( s * s + 1.0 );
                    e[i] = g * r;
                    c = 1.0 / r;
                    s = s * c;
                }
                g = d[i] - p;
                r = ( d[i-1] - g ) * s + 2.0 * c * b;
                p = s * r;
                d[i] = g + p;
                g = c * r - b;
                f = z[i];
                z[i] = s * z[i-1] + c * f;
                z[i-1] = c * z[i-1] - s * f;
            }
            d[l-1] = d[l-1] - p;
            e[l-1] = g;
            e[m-1] = 0.0;
        }
    }

    //  Sorting.
    for (int ii = 2; ii <= m; ii++ ) {
        int i = ii - 1;
        int k = i;
        fpv p = d[i-1];

        for (int j = ii; j <= n; j++ ) {
            if ( d[j-1] < p ){
                k = j;
                p = d[j-1];
            }
        }

        if ( k != i ) {
            d[k-1] = d[i-1];
            d[i-1] = p;
            p = z[i-1];
            z[i-1] = z[k-1];
            z[k-1] = p;
        }
    }
}
//****************************************************************************80


//****************************************************************************80
//
//  Purpose:
//
//    SGQF computes knots and weights of a Gauss Quadrature formula.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with simple knots from the Jacobi matrix and the zero-th
//    moment of the weight function, using the Golub-Welsch technique.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//    Input, double AJ[NT], the diagonal of the Jacobi matrix.
//    Input/output, double BJ[NT], the subdiagonal of the Jacobi
//        matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
//    Input, double ZEMU, the zero-th moment of the weight function.
//    Output, double T[NT], the knots.
//    Output, double WTS[NT], the weights.
//
//****************************************************************************80
template<typename fpv>
void sgqf ( int nt, fpv* aj, fpv* bj, fpv zemu, fpv* t, fpv* wts) {
    //  Exit if the zero-th moment is not positive.
    if ( zemu <= 0.0 ) crash("SGQF - Fatal error! ZEMU <= 0.");

    //  Set up vectors for IMTQLX.
    for (int i = 0; i < nt; i++ ) t[i] = aj[i];
    wts[0] = sqrt ( zemu );
    for (int i = 1; i < nt; i++ ) wts[i] = 0.0;

    //  Diagonalize the Jacobi matrix.
    imtqlx ( nt, t, bj, wts );

    for (int i = 0; i < nt; i++ ) wts[i] *= wts[i];
}
//****************************************************************************80


//****************************************************************************80
//
//  Purpose:
//
//    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with a classical weight function with default values for A and B,
//    and only simple knots.
//
//    There are no moments checks and no printing is done.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, double ALPHA, the value of Alpha, if needed.
//    Input, double BETA, the value of Beta, if needed.
//    Output, double T[NT], the knots.
//    Output, double WTS[NT], the weights.
//
//****************************************************************************80

template<typename fpv>
void cdgqf (int nt, int kind, fpv alpha, fpv beta, fpv* t, fpv* wts) {
    parchk ( kind, 2 * nt, alpha, beta );

    //  Get the Jacobi matrix and zero-th moment.
    fpv* aj = new fpv[nt];
    fpv* bj = new fpv[nt];
    fpv zemu = class_matrix ( kind, nt, alpha, beta, aj, bj );

    //  Compute the knots and weights.
    sgqf ( nt, aj, bj, zemu, t, wts );

    delete [] aj;
    delete [] bj;
}
//****************************************************************************80
