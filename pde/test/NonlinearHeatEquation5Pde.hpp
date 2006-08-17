#ifndef _NONLINEARHEATEQUATION5PDE_HPP_
#define _NONLINEARHEATEQUATION5PDE_HPP_

#include "AbstractNonlinearEllipticPde.hpp"


template <int SPACE_DIM>
class NonlinearHeatEquation5Pde : public AbstractNonlinearEllipticPde<SPACE_DIM>
{
public:

    double ComputeLinearSourceTerm(Point<SPACE_DIM> )
    {
        return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double u)
    {
        return (-3*u*u*u);
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> , double u)
    {
        return identity_matrix<double>(SPACE_DIM) * (u*u);
    }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(Point<SPACE_DIM> , double u)
    {
        return identity_matrix<double>(SPACE_DIM) * 2.0*u;
    }
    
    double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> , double u)
    {
        return (-9*u*u);//(-(-4*exp(-x[0])+4*u));
    }
};

#endif //_NONLINEARHEATEQUATION5PDE_HPP_
