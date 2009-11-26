/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef NHSMODELWITHIMPLICITSOLVER_HPP_
#define NHSMODELWITHIMPLICITSOLVER_HPP_

#include "NhsContractionModel.hpp"


/**
 *  NHS system with built-in backward solver. Jon Whiteley's method, which breaks down
 *  the multivariable implicit solve into a sequence of 1d implicit solves.
 *
 *  The full backward Euler method on the NHS system compute the solution w^{n+1} \in R^5
 *  by:
 *   w^{n+1}  - dt g(w^{n+1})  =  w^n
 *  where the ODE system is dw/dt = g(w).
 * 
 *  This would involve solving a 5d nonlinear system at each timestep. However, only g1 and g2
 *  (ie dCatrop/dt and dz/dt) are nonlinear, g3-g5 are linear and uncoupled. Therefore 
 *  the backward euler solutions of Q1,Q2,Q3 can be computed immediately, leaving a 2D nonlinear
 *  system.   
 *  
 *  This system could be solved with newton's method as normal, here an alternative method  
 *  in which it is broken down into a sequence of 1D problems is used. 
 * 
 *  The idea is to solve a 1d nonlinear problem for the active tension, where the problem
 *  f(T_a)=T_a, where T_a is the current active tension, and f(T_a) is the active tension
 *  obtained by assuming this active tension and solving the ODEs as follows:
 *
 *  Given T_a, so compute solution for Ca_trop using backward euler (1d linear problem)
 *  Using this Ca_trop, compute solution for z using backward euler (1d nonlinear problem)
 *  Solve for Qi implicitly using backward euler (three 1d linear problems)
 *  Compute T0, Q and then T_a
 *
 *  Reference: J.P. Whiteley, M.J. Bishop, D.J. Gavaghan "Soft Tissue Modelling of Cardiac
 *  Fibres for Use in Coupled Mechano-Electric Simulations"
 */
class NhsModelWithImplicitSolver : public NhsContractionModel
{
private:
    /**
     *  Tolerance for solving nonlinear system which require newton's method
     *  NOTE: 1e-6 doesn't give graphically close results when solving a full
     *  problem and comparing implicit with explicit
     */
    const static double mTolerance = 1e-10;

    /** Timestep for the ODEs solving */
    double mDt;

    /** See SetUseImplicitExplicitSolveForZ() */
    bool mUseImplicitExplicitSolveForZ;


    /** Current state variables to be used in the next timestep */
    std::vector<double> mCurrentStateVars;

    /** Initial guess for the active tension */
    double mActiveTensionInitialGuess;
    /** The solution for the active tension after a RunDoNotUpdate() has been called */
    double mActiveTensionSolution;

    /**
     *  The main Newton solve
     *
     *  Solve a 1d nonlinear system, f(T_a)=T_a, where T_a is the current active tension
     *  and f(T_a) is the active tension obtained by assuming this active tension and solving
     *  the ODEs as follows
     *
     *  Given T_a, so compute solution for Ca_trop using backward euler (1d linear problem)
     *  Using this Ca_trop, compute solution for z using backward euler (1d nonlinear problem)
     *  Solve for Qi implicitly using backward euler (three 1d linear problems)
     *  Compute T0, Q and then T_a
     */
    void ImplicitSolveForActiveTension();
    
    /** The residual function for the main Newton solve. See ImplicitSolveForActiveTension()
     * 
     * @param activeTensionGuess
     * @param Q = Q1+Q2+Q3 which can be solved for first
     */
    double CalcActiveTensionResidual(double activeTensionGuess, double Q);

    /**
     *  Assume the active tension is known and solve for the Ca_trop at the next time
     *  implicitly using backward euler. This can be done directly as the rhs is linear
     *  in Ca_trop
     * 
     * @param newActiveTension
     */
    double ImplicitSolveForCaTrop(double newActiveTension);

    /**
     *  Assume the Ca_trop is known and solve for the z at the next time
     *  implicitly using backward euler. Uses Newton's method.
     * 
     * @param newCaTrop
     */
    double ImplicitSolveForZ(double newCaTrop);
    /**
     * Residual for solving z implicitly. See ImplicitSolveForZ().
     * 
     * @param z
     * @param newCaTrop
     */
    double CalcZResidual(double z, double newCaTrop);
    /** Solve for z semi implicitly. See UseImplicitExplicitSolveForZ().
     * 
     * @param newCaTrop
     */
    double ImplicitExplicitSolveForZ(double newCaTrop);

    /**
     *  Solve for Q1,Q2,Q3 (and therefore Q) implicitly using backward euler.
     *  These can be done directly as the rhs is linear in Qi
     */
    double ImplicitSolveForQ();


public :
    /**
     *  Constructor
     */
    NhsModelWithImplicitSolver();

    /**
     *  Set a current active tension guess. Generally not needed as the current
     *  active tension is used if this isn't called.
     * 
     *  @param activeTensionInitialGuess
     */
    void SetActiveTensionInitialGuess(double activeTensionInitialGuess);

    /**
     *  Solves for the new state variables at the given end time using the implicit
     *  method. Note that the internal state variables are not altered, the solution
     *  is saved instead. Call UpdateStateVariables() to update, and
     *  GetNextActiveTension() to get the solved active tension
     *
     *  The state variables are not updated because this solve will be called as part
     *  of the newton iteration (ie guess stretch, see what the new active tension is)
     *  in a fully implicit method.
     * 
     *  Note: overloaded from the method in AbstractOdeBasedContractionModel, which
     *  just does a simple Euler solve
     * 
     *  @param startTime
     *  @param endTime
     *  @param timestep
     */
    void RunDoNotUpdate(double startTime, double endTime, double timestep);

    /**
     *  Solve for z semi-implicitly instead of fully implicitly. If we assume we know
     *  Ca_trop solving for z is a 1d nonlinear problem. Call this to treat the problem
     *  implicitly in the linear terms on the rhs of dzdt (the (1-z) and (z) terms), and
     *  explicitly in the nonlinear term (the z^nr/(z^nr + K^nr) term. This means the
     *  problem can be solved directly and no Newton iterations are needed.
     * 
     * @param useImplicitExplicitSolveForZ
     */
    void UseImplicitExplicitSolveForZ(bool useImplicitExplicitSolveForZ = true);

    /**
     *  Get the active tension corresponding to the stored state variables computed
     *  from the last RunDoNotUpdate(), ie the active tension at the next time.
     *  Note that calling GetActiveTension() on the base class will use the internal
     *  state variables and return the active tension at the last time, if
     *  RunDoNotUpdate() has been called but UpdateStateVariables() has not
     */
    double GetNextActiveTension();
};

#endif /*NHSMODELWITHIMPLICITSOLVER_HPP_*/
