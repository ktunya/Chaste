/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef UBLASCUSTOMFUNCTIONS_HPP_
#define UBLASCUSTOMFUNCTIONS_HPP_

/**
 * @file
 * A collection of useful functions extending the functionality of the
 * Boost Ublas library.
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include "petsc.h"
#include "petscblaslapack.h"
//Promote universal LAPACK name if it's an old version of PETSc
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
#define LAPACKgeev_ LAgeev_
#endif

#include "Exception.hpp"
#include "MathsCustomFunctions.hpp"
#include "PetscTools.hpp"

using namespace boost::numeric::ublas;

// COMMON DETERMINANTS - SQUARE

/**
 * 1x1 Determinant.
 * Get the determinant of a ublas matrix.
 *
 * @param rM The matrix of which to find the determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM)
{
    using namespace boost::numeric::ublas;

    return rM(0,0);
}

/**
 * 2x2 Determinant.
 * Get the determinant of a ublas matrix.
 *
 * @param rM The matrix of which to find the determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T,2,2>& rM)
{
    using namespace boost::numeric::ublas;

    return rM(0,0)*rM(1,1) - rM(1,0)*rM(0,1);
}

/**
 * 3x3 Determinant.
 * Get the determinant of a ublas matrix.
 *
 * @param rM The matrix of which to find the determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM)
{
    using namespace boost::numeric::ublas;

    return    rM(0,0) * (rM(1,1)*rM(2,2) - rM(1,2)*rM(2,1))
            - rM(0,1) * (rM(1,0)*rM(2,2) - rM(1,2)*rM(2,0))
            + rM(0,2) * (rM(1,0)*rM(2,1) - rM(1,1)*rM(2,0));
}

// COMMON GENERALIZED DETERMINANTS - NOT SQUARE

/**
 * 3x2 (Generalized determinant).
 * Calculate the generalized determinant of a 3x2 matrix.
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM)
{
    using namespace boost::numeric::ublas;
    c_matrix<T,2,2> product = prod(trans(rM), rM);
    return std::sqrt(Determinant(product));
}

/**
 * 3x1 (Generalized determinant).
 * Calculate the generalized determinant of a 3x1 matrix.
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM)
{
    using namespace boost::numeric::ublas;
    return std::sqrt(rM(0,0)*rM(0,0) + rM(1,0)*rM(1,0) + rM(2,0)*rM(2,0));
}

/**
 * 2x1 (Generalized determinant).
 * Calculate the generalized determinant of a 2x1 matrix.
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM)
{
    using namespace boost::numeric::ublas;
    return   std::sqrt(rM(0,0) * rM(0,0) + rM(1,0) * rM(1,0));
}

/**
 * 3x0 (Generalized determinant) - not implement, but needed by some compilers.
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 0>& rM)
{
    NEVER_REACHED;
}

/**
 * 2x0 (Generalized determinant) - not implement, but needed by some compilers.
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 2, 0>& rM)
{
    NEVER_REACHED;
}

/**
 * 1x0 (Generalized determinant) - not implement, but needed by some compilers.
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 1, 0>& rM)
{
    NEVER_REACHED;
}

// COMMON SUBDETERMINANTS - SQUARE

/**
 * 1x1 SubDeterminant.
 * Return the determinant of a submatrix after removing a particular row and column
 * For a 1x1 matrix this should always remove the only row and column (0,0).
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow == 0);
    assert(misscol == 0);
    return 1.0;
}

/**
 * 2x2 SubDeterminant.
 * Return the determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 2>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 2);
    assert(misscol < 2);

    unsigned row = (missrow==1) ? 0 : 1;
    unsigned col = (misscol==1) ? 0 : 1;
    return rM(row,col);
}

/**
 * SubDeterminant 3x3.
 * Determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    assert(misscol < 3);

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

// COMMON SUBDETERMINANTS - NOT SQUARE

/**
 * SubDeterminant 3x2.
 * Determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    //assert(misscol < 2); //Don't assert this since it is used

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

/**
 * SubDeterminant 3x1.
 * Determinant of a submatrix after removing a particular row and column.
 * @param rM The matrix of which to find the subdeterminant.
 *
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    assert(misscol < 1);

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

/**
 * SubDeterminant 2x1.
 * Determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 2);
    assert(misscol < 1);

    unsigned row = (missrow==1) ? 0 : 1;
    unsigned col = (misscol==1) ? 0 : 1;
    return rM(row,col);
}

#if defined(__xlC__)
/* IBM compiler doesn't support zero-sized arrays*/
#else //#if defined(__xlC__)
/**
 * SubDeterminant 3x0 - Not implemented, but needed by some compilers for recursive template calls.
 * Determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}

/**
 * SubDeterminant 2x0 - Not implemented, but needed by some compilers for recursive template calls.
 * Determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}

/**
 * SubDeterminant 1x0 - Not implemented, but needed by some compilers for recursive template calls.
 * Determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 1, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}
#endif //#if defined(__xlC__)

// COMMON INVERSES - SQUARE

/**
 * 1x1 Inverse.
 * Get the inverse of a ublas matrix.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 1> Inverse(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T,1,1> inverse;
    T det = Determinant(rM);
    assert(fabs(det) > DBL_EPSILON); // else it is a singular matrix
    inverse(0,0) =  1.0/det;
    return inverse;
}

/**
 * 2x2 Inverse.
 * Get the inverse of a ublas matrix.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 2, 2> Inverse(const boost::numeric::ublas::c_matrix<T, 2, 2>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 2, 2> inverse;
    T det = Determinant(rM);

    assert( fabs(det) > DBL_EPSILON ); // else it is a singular matrix
    inverse(0,0)  =  rM(1,1)/det;
    inverse(0,1)  = -rM(0,1)/det;
    inverse(1,0)  = -rM(1,0)/det;
    inverse(1,1)  =  rM(0,0)/det;
    return inverse;
}

/**
 * 3x3 Inverse.
 * Get the inverse of a ublas matrix.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 3, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 3, 3> inverse;
    T det = Determinant(rM);
    assert(fabs(det) > DBL_EPSILON); // else it is a singular matrix

    inverse(0,0) =  (rM(1,1)*rM(2,2) - rM(1,2)*rM(2,1))/det;
    inverse(1,0) = -(rM(1,0)*rM(2,2) - rM(1,2)*rM(2,0))/det;
    inverse(2,0) =  (rM(1,0)*rM(2,1) - rM(1,1)*rM(2,0))/det;
    inverse(0,1) = -(rM(0,1)*rM(2,2) - rM(0,2)*rM(2,1))/det;
    inverse(1,1) =  (rM(0,0)*rM(2,2) - rM(0,2)*rM(2,0))/det;
    inverse(2,1) = -(rM(0,0)*rM(2,1) - rM(0,1)*rM(2,0))/det;
    inverse(0,2) =  (rM(0,1)*rM(1,2) - rM(0,2)*rM(1,1))/det;
    inverse(1,2) = -(rM(0,0)*rM(1,2) - rM(0,2)*rM(1,0))/det;
    inverse(2,2) =  (rM(0,0)*rM(1,1) - rM(0,1)*rM(1,0))/det;

    return inverse;
}

// COMMON PSEUDO-INVERSES - NOT SQUARE

/**
 * 2x3 pseudo-inverse of a  matrix.
 * The pseudo-inverse is given by pinv(T) = (T'T)^(-1)*T'.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 2, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 2, 3> inverse;

    //
    // calculate (T'T)^-1, where T'T = (a b)
    //                                 (c d)

    T a = rM(0,0)*rM(0,0) + rM(1,0)*rM(1,0) + rM(2,0)*rM(2,0);
    T b = rM(0,0)*rM(0,1) + rM(1,0)*rM(1,1) + rM(2,0)*rM(2,1);
    T c = b;
    T d = rM(0,1)*rM(0,1) + rM(1,1)*rM(1,1) + rM(2,1)*rM(2,1);

    T det = a*d - b*c;

    T a_inv =  d/det;
    T b_inv = -b/det;
    T c_inv = -c/det;
    T d_inv =  a/det;

    inverse(0,0) = a_inv*rM(0,0) + b_inv*rM(0,1);
    inverse(1,0) = c_inv*rM(0,0) + d_inv*rM(0,1);
    inverse(0,1) = a_inv*rM(1,0) + b_inv*rM(1,1);
    inverse(1,1) = c_inv*rM(1,0) + d_inv*rM(1,1);
    inverse(0,2) = a_inv*rM(2,0) + b_inv*rM(2,1);
    inverse(1,2) = c_inv*rM(2,0) + d_inv*rM(2,1);

    return inverse;
}

/**
 * 2x1 pseudo-inverse of a  matrix.
 * The pseudo-inverse is given by pinv(T) = (T'T)^(-1)*T'.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 2> Inverse(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 1, 2> inverse;
    T det = Determinant(rM);

    inverse(0,0) = rM(0,0)/det/det;
    inverse(0,1) = rM(1,0)/det/det;

    return inverse;
}

/**
 * 3x1 pseudo-inverse of a  matrix.
 * The pseudo-inverse is given by pinv(T) = (T'T)^(-1)*T'.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 1, 3> inverse;
    T det = Determinant(rM);

    inverse(0,0) = rM(0,0)/det/det;
    inverse(0,1) = rM(1,0)/det/det;
    inverse(0,2) = rM(2,0)/det/det;

    return inverse;
}

// COMMON MATRIX TRACES

/**
 * 1x1 matrix trace (sum of diagonal elements).
 *
 * @param rM The matrix of which to find the trace.
 * @return The trace
 */
template<class T>
T Trace(const c_matrix<T, 1, 1>& rM)
{
    return rM(0,0);
}

/**
 * 2x2 matrix trace (sum of diagonal elements).
 *
 * @param rM The matrix of which to find the trace.
 * @return The trace
 */
template<class T>
T Trace(const c_matrix<T, 2, 2>& rM)
{
    return rM(0,0) + rM(1,1);
}

/**
 * 3x3 matrix trace (sum of diagonal elements).
 *
 * @param rM The matrix of which to find the trace.
 * @return The trace
 */
template<class T>
T Trace(const c_matrix<T, 3, 3>& rM)
{
    return rM(0,0) + rM(1,1) + rM(2,2);
}

/**
 * 4x4 matrix trace (sum of diagonal elements).
 *
 * @param rM The matrix of which to find the trace.
 * @return The trace
 */
template<class T>
T Trace(const c_matrix<T, 4, 4>& rM)
{
    return rM(0,0) + rM(1,1) + rM(2,2) + rM(3,3);
}

// OTHER MATRIX FUNCTIONS (INVARIANTS, EIGENVECTORS)

/**
 * 3x3 second invariant.
 * @note Implementation only correct for
 * a SYMMETRIC matrix though. It is up to the user to check the
 * input matrix is symmetric.
 *
 * @param rM The matrix
 */
template<class T>
T SecondInvariant(const c_matrix<T, 3, 3>& rM)
{
    return    rM(0,0)*rM(1,1) + rM(1,1)*rM(2,2) + rM(2,2)*rM(0,0)
            - rM(1,0)*rM(1,0) - rM(2,1)*rM(2,1) - rM(2,0)*rM(2,0);
}

/**
 * 2x2 second invariant.
 * Second invariant of a 2d matrix, i.e. the determinant. This function
 * is mainly here just so that the same code can be used in 2d and 3d.
 *
 * @param rM The matrix
 */
template<class T>
T SecondInvariant(const c_matrix<T, 2, 2>& rM)
{
    return Determinant(rM);
}

/**
 * Use LAPACK functionality to find the eigenvector corresponding
 * real eigenvalue which is smallest in magnitude.
 * Caveat: if there are zero eigenvalues they are ignored.
 * It's the smallest magnitude non-zero real eigenvalue which is used.
 *
 * @param rA 3x3 matrix is question
 * @return 3-vector corresponding to right-eigenvector in question
 */
c_vector<double,3> CalculateEigenvectorForSmallestNonzeroEigenvalue(c_matrix<double, 3, 3>& rA);

//COMMON VECTOR FUNCTIONS

/**
 * This is a cross-product aka vector-product, only implemented for 3-vectors.
 *
 * @param rA first vector
 * @param rB second vector
 * @return rA x rB
 */
template<class T>
c_vector<T, 3> VectorProduct(const c_vector<T, 3>& rA, const c_vector<T, 3>& rB)
{

    c_vector<T, 3> result;

    double x1 = rA(0);
    double y1 = rA(1);
    double z1 = rA(2);
    double x2 = rB(0);
    double y2 = rB(1);
    double z2 = rB(2);

    result(0) = y1*z2 - z1*y2;
    result(1) = z1*x2 - x1*z2;
    result(2) = x1*y2 - y1*x2;

    return result;
}

/**
 * Convenience function for quickly creating test vectors (1D).
 *
 * @param x entry in vector
 * @returns vector=(x)
 */
c_vector<double, 1> Create_c_vector(double x);

/**
 * Convenience function for quickly creating test vectors (2D).
 *
 * @param x entry in vector
 * @param y entry in vector
 * @returns vector=(x,y)
 */
c_vector<double, 2> Create_c_vector(double x, double y);

/**
 * Convenience function for quickly creating test vectors (3D).
 *
 * @param x entry in vector
 * @param y entry in vector
 * @param z entry in vector
 * @returns vector=(x,y,z)
 */
c_vector<double, 3> Create_c_vector(double x, double y, double z);

#endif /*UBLASCUSTOMFUNCTIONS_HPP_*/
