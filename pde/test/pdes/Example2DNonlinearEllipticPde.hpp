/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef _EXAMPLE2DNONLINEARELLIPTICPDE_HPP_
#define _EXAMPLE2DNONLINEARELLIPTICPDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

/**
 * A relatively simple PDE for testing the nonlinear elliptic assembler in 2D.
 */

class Example2DNonlinearEllipticPde:public AbstractNonlinearEllipticPde<2>
{

public:
    double ComputeLinearSourceTerm(const ChastePoint<2>& p)
    {
        double x = p[0];
        double y = p[1];
        return -( 4 + 8*x*x + 8*y*y );
    }

    double ComputeNonlinearSourceTerm(const ChastePoint<2>& , double )
    {
        return 0.0;
    }

    c_matrix<double, 2, 2> ComputeDiffusionTerm(const ChastePoint<2>& , double u)
    {
        return identity_matrix<double>(2)*u;
    }

    c_matrix<double, 2, 2> ComputeDiffusionTermPrime(const ChastePoint<2>& , double )
    {
        return identity_matrix<double>(2);
    }

    double ComputeNonlinearSourceTermPrime(const ChastePoint<2>& , double )
    {
        return 0.0;
    }

    virtual ~Example2DNonlinearEllipticPde()
    {}
};

#endif //_EXAMPLE2DNONLINEARELLIPTICPDE_HPP_
