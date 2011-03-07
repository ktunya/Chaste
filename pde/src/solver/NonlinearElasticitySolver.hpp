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
#ifndef NONLINEARELASTICITYSOLVER_HPP_
#define NONLINEARELASTICITYSOLVER_HPP_

/*
 * NOTE ON COMPILATION ERRORS:
 *
 * This file won't compile with Intel icpc version 9.1.039, with error message:
 * "Terminate with:
  (0): internal error: backend signals"
 *
 * Try recompiling with icpc version 10.0.025.
 */


#include "AbstractNonlinearElasticitySolver.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"

/**
 *  Finite elasticity solver. Solves static *incompressible* nonlinear elasticity
 *  problems with arbitrary (incompressible) material laws and a body force.
 *
 *  Uses quadratic-linear bases (for displacement and pressure), and is therefore
 *  outside other assembler or solver hierarchy.
 */
template<size_t DIM>
class NonlinearElasticitySolver : public AbstractNonlinearElasticitySolver<INCOMPRESSIBLE,DIM>
{
    friend class TestNonlinearElasticitySolver;
    friend class TestNonlinearElasticityAdjointSolver;
    friend class AdaptiveNonlinearElasticityProblem;

protected:

    /** Number of vertices per element */
    static const size_t NUM_VERTICES_PER_ELEMENT = DIM+1;
    /** Number of nodes per element */
    static const size_t NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic
    /** Stencil size - number of unknowns per element (DIM*NUM_NODES_PER_ELEMENT displacement unknowns,
     *  NUM_VERTICES_PER_ELEMENT pressure unknowns */
    static const size_t STENCIL_SIZE = DIM*NUM_NODES_PER_ELEMENT + NUM_VERTICES_PER_ELEMENT;
    /** Number of nodes per boundary element */
    static const size_t NUM_NODES_PER_BOUNDARY_ELEMENT = DIM*(DIM+1)/2;
    /** Boundary stencil size */
    static const size_t BOUNDARY_STENCIL_SIZE = DIM*NUM_NODES_PER_BOUNDARY_ELEMENT + DIM;



    /**
     *  The material laws for each element. This will either be of size 1 (same material law for all elements,
     *  ie homogeneous), or size num_elem.
     */
    std::vector<AbstractIncompressibleMaterialLaw<DIM>*> mMaterialLaws;

    /**
     *  The solution pressures. mPressures[i] = pressure at node i (ie vertex i).
     */
    std::vector<double> mPressures;


    /**
     * Assemble residual or jacobian on an element, using the current solution
     * stored in mCurrrentSolution. The ordering assumed is (in 2d)
     * rBElem = [u0 v0 u1 v1 .. u5 v5 p0 p1 p2].
     *
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     * @param rAElemPrecond The element's contribution to the matrix passed to PetSC
     *     in creating a preconditioner
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    virtual void AssembleOnElement(Element<DIM, DIM>& rElement,
                                   c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                                   c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElemPrecond,
                                   c_vector<double, STENCIL_SIZE>& rBElem,
                                   bool assembleResidual,
                                   bool assembleJacobian);

    /**
     * Compute the term from the surface integral of s*phi, where s is
     * a specified non-zero surface traction (ie Neumann boundary condition)
     * to be added to the Rhs vector.
     *
     * @param rBoundaryElement
     * @param rAelem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     * @param rBelem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param rTraction
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    virtual void AssembleOnBoundaryElement(BoundaryElement<DIM-1, DIM>& rBoundaryElement,
                                           c_matrix<double, BOUNDARY_STENCIL_SIZE, BOUNDARY_STENCIL_SIZE>& rAelem,
                                           c_vector<double, BOUNDARY_STENCIL_SIZE>& rBelem,
                                           c_vector<double, DIM>& rTraction,
                                           bool assembleResidual,
                                           bool assembleJacobian);

    /**
     *  Set up the current guess to be the solution given no displacement.
     *  The current solution (in 2d) is order as
     *  [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *  (where there are N total nodes and M vertices)
     *  so the initial guess is
     *  [0 0 0 0 ... 0 0 p1 p2 .. pM]
     *  where p_i are such that T is zero (depends on material law).
     *
     *  In a homogeneous problem, all p_i are the same.
     *  In a heterogeneous problem, p for a given vertex is the
     *  zero-strain-pressure for ONE of the elements containing that
     *  vertex (which element containing the vertex is reached LAST). In
     *  this case the initial guess will be close but not exactly the
     *  solution given zero body force.
     */
    void FormInitialGuess();

    /**
     * Assemble the residual vector (using the current solution stored
     * in mCurrentSolution, output going to mpLinearSystem->rGetRhsVector),
     * or Jacobian matrix (using the current solution stored in
     * mCurrentSolution, output going to mpLinearSystem->rGetLhsMatrix).
     *
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    void AssembleSystem(bool assembleResidual, bool assembleJacobian);

    
    /**
     *  Simple (one-line function which just calls ComputeStressAndStressDerivative() on the material law, using C,
     *  inv(C), and p as the input and with rT and rDTdE as the output. Overloaded by other assemblers (eg cardiac 
     *  mechanics) which need to add extra terms to the stress.
     * 
     *  @param pMaterialLaw The material law for this element
     *  @param rC The Lagrangian deformation tensor (F^T F)
     *  @param rInvC The inverse of C. Should be computed by the user.
     *  @param pressure The current pressure
     *  @param elementIndex Index of the current element 
     *  @param currentQuadPointGlobalIndex The index (assuming an outer loop over elements and an inner 
     *    loop over quadrature points), of the current quadrature point.
     *  @param rT The stress will be returned in this parameter
     *  @param rDTdE the stress derivative will be returned in this parameter, assuming
     *    the final parameter is true
     *  @param computeDTdE A boolean flag saying whether the stress derivative is
     *    required or not.
     */ 
    virtual void ComputeStressAndStressDerivative(AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                                  c_matrix<double,DIM,DIM>& rC, 
                                                  c_matrix<double,DIM,DIM>& rInvC,
                                                  double pressure,
                                                  unsigned elementIndex,
                                                  unsigned currentQuadPointGlobalIndex,
                                                  c_matrix<double,DIM,DIM>& rT,
                                                  FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                                  bool computeDTdE);

public:

    /**
     * Constructor for homogeneous problems.
     *
     * @param pQuadMesh The quadratic mesh to solve on
     * @param pMaterialLaw A single material law to use on all elements
     * @param bodyForce The body force if constant. (If not constant, pass in a zero vector and call SetFunctionalBodyForce()
     * @param density The density (assumed constant)
     * @param outputDirectory The output directory
     * @param fixedNodes Which nodes are fixed in space (the displacement is assumed to be zero unless the next parameter is given
     * @param pFixedNodeLocations Optional new locations of the fixed nodes.
     */
    NonlinearElasticitySolver(QuadraticMesh<DIM>* pQuadMesh,
                              AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                              c_vector<double,DIM> bodyForce,
                              double density,
                              std::string outputDirectory,
                              std::vector<unsigned>& fixedNodes,
                              std::vector<c_vector<double,DIM> >* pFixedNodeLocations = NULL);

    /**
     * Variant constructor taking a vector of material laws for heterogeneous problems
     *
     * @param pQuadMesh The quadratic mesh to solve on
     * @param rMaterialLaws Vector of material laws for each element
     * @param bodyForce The body force if constant. (If not constant, pass in a zero vector and call SetFunctionalBodyForce()
     * @param density The density (assumed constant)
     * @param outputDirectory The output directory
     * @param fixedNodes Which nodes are fixed in space (the displacement is assumed to be zero unless the next parameter is given
     * @param pFixedNodeLocations Optional new locations of the fixed nodes.
     */
    NonlinearElasticitySolver(QuadraticMesh<DIM>* pQuadMesh,
                              std::vector<AbstractIncompressibleMaterialLaw<DIM>*>& rMaterialLaws,
                              c_vector<double,DIM> bodyForce,
                              double density,
                              std::string outputDirectory,
                              std::vector<unsigned>& fixedNodes,
                              std::vector<c_vector<double,DIM> >* pFixedNodeLocations = NULL);

    /** Destructor. */
    ~NonlinearElasticitySolver();

    /**
     * Get pressures for each vertex
     */
    std::vector<double>& rGetPressures();

};

#endif /*NONLINEARELASTICITYSOLVER_HPP_*/
