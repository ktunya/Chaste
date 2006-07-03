#ifndef _BIDOMAINDG0ASSEMBLER_HPP_
#define _BIDOMAINDG0ASSEMBLER_HPP_

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "Point.hpp"
#include "Element.hpp"
#include "LinearSystem.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "AbstractLinearSolver.hpp"
#include "BidomainPde.hpp"
#include "AbstractBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "AbstractAssembler.hpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class BidomainDg0Assembler : public AbstractAssembler<ELEMENT_DIM,SPACE_DIM>
{
private:
    // pde
    BidomainPde<SPACE_DIM>* mpBidomainPde;
    
    // linear system and solver
    LinearSystem* mpAssembledLinearSystem;
    AbstractLinearSolver* mpSolver;
    bool mMatrixIsAssembled;
       
    // mesh
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* mpMesh;    
        
    // times
    double mTstart;
    double mTend;
    double mDt;
    bool   mTimesSet;

    // initial condition
    Vec  mInitialCondition;
    bool mInitialConditionSet;
    
    
    void AssembleOnElement(const Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                           c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2>& rAElem,
                           c_vector<double, 2*ELEMENT_DIM+2>& rBElem)
    {
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpQuadRule);
        AbstractBasisFunction<ELEMENT_DIM> &rBasisFunction =
            *(AbstractAssembler<ELEMENT_DIM,SPACE_DIM>::mpBasisFunction);
        
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *inverseJacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if (!mMatrixIsAssembled)
        {
            inverseJacobian = rElement.GetInverseJacobian();
            rAElem.clear(); 
        }
        rBElem.clear();

                        
        const int num_elem_nodes = rElement.GetNumNodes();
                
        // loop over guass points
        for (int quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            Point<ELEMENT_DIM> quad_point = quad_rule.GetQuadPoint(quad_index);

            c_vector<double, ELEMENT_DIM+1> basis_func = rBasisFunction.ComputeBasisFunctions(quad_point);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>  grad_basis;
            
            if (!mMatrixIsAssembled)
            {
                grad_basis = rBasisFunction.ComputeTransformedBasisFunctionDerivatives
                                                (quad_point, *inverseJacobian);
            }
            Point<SPACE_DIM> x(0,0,0);

            // Vm is the trans-membrane voltage interpolated onto the gauss point
            // I_ionic is the ionic current (per unit AREA) interpolated onto gauss point
            // I_intra_stim, I_extra_stim are the stimuli (per unit VOLUME) interpolated 
            //  onto gauss point
            double Vm = 0;
            double I_ionic = 0;
            double I_intra_stim = 0;
            double I_extra_stim = 0;

            // interpolate x, Vm, and currents
            for (int i=0; i<num_elem_nodes; i++)
            {
                const Node<SPACE_DIM> *node = rElement.GetNode(i);
                const Point<SPACE_DIM> node_loc = node->rGetPoint();
                for (int j=0; j<SPACE_DIM; j++)
                {
                    x.SetCoordinate(j, x[j] + basis_func(i)*node_loc[j]);
                }
                
                int node_global_index = rElement.GetNodeGlobalIndex(i);
                
                Vm           += basis_func(i)*mpBidomainPde->GetInputCacheMember( 2*node_global_index );
                I_ionic      += basis_func(i)*mpBidomainPde->GetIionicCacheReplicated()[ node_global_index ];
                I_intra_stim += basis_func(i)*mpBidomainPde->GetIntracellularStimulusCacheReplicated()[ node_global_index ];
                I_extra_stim += basis_func(i)*mpBidomainPde->GetExtracellularStimulusCacheReplicated()[ node_global_index ];
            }

            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);

            // get bidomain parameters
            double Am = mpBidomainPde->GetSurfaceAreaToVolumeRatio();
            double Cm = mpBidomainPde->GetCapacitance();
            
            c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_i = mpBidomainPde->GetIntracellularConductivityTensor();
            c_matrix<double, SPACE_DIM, SPACE_DIM> sigma_e = mpBidomainPde->GetExtracellularConductivityTensor();
            

            // assemble element stiffness matrix
            if (!mMatrixIsAssembled)
            {
                c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp = prod(sigma_i, grad_basis);
                c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_i_grad_phi = 
                    prod(trans(grad_basis), temp);
                                                     
                c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> basis_outer_prod =
                    outer_prod(basis_func, basis_func);
                
                c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> temp2 = prod(sigma_e, grad_basis);
                c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> grad_phi_sigma_e_grad_phi = 
                    prod(trans(grad_basis), temp2);
                                           
                // Components of the element stiffness matrix are:   
                // (0,0) block:            ACV/dt + (Di grad_basis_col)dot(grad_basis_row)
                // (0,1) and (1,0) blocks: (Di grad_basis_col)dot(grad_basis_row)
                // (1,1) block:           ( ((Di+De)grad_basis_col )dot(grad_basis_row) 
                
                /* old version:
                rAElem(2*row,  2*col)   += wJ*( (Am*Cm/mDt)*basis_func(col)*basis_func(row)  +  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col )) );
                rAElem(2*row+1,2*col)   += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))   );
                rAElem(2*row,  2*col+1) += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))  );
                rAElem(2*row+1,2*col+1) += wJ*(  inner_prod( grad_basis_row, prod( sigma_i, grad_basis_col ))   +   inner_prod( grad_basis_row, prod( sigma_e, grad_basis_col )));
                */
                
                // a matrix slice is a submatrix
                
                // even rows, even columns
                matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
                    rAElem_slice00(rAElem, slice (0, 2, num_elem_nodes), slice (0, 2, num_elem_nodes));
                rAElem_slice00 += wJ*( (Am*Cm/mDt)*basis_outer_prod + grad_phi_sigma_i_grad_phi );
                
                // odd rows, even columns
                matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
                    rAElem_slice10(rAElem, slice (1, 2, num_elem_nodes), slice (0, 2, num_elem_nodes));
                rAElem_slice10 += wJ * grad_phi_sigma_i_grad_phi;

                // even rows, odd columns
                matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
                    rAElem_slice01(rAElem, slice (0, 2, num_elem_nodes), slice (1, 2, num_elem_nodes));
                rAElem_slice01 += wJ * grad_phi_sigma_i_grad_phi;
                
                // odd rows, odd columns
                matrix_slice<c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> >
                    rAElem_slice11(rAElem, slice (1, 2, num_elem_nodes), slice (1, 2, num_elem_nodes));
                rAElem_slice11 += wJ*( grad_phi_sigma_i_grad_phi + grad_phi_sigma_e_grad_phi);
            }
            
            // assemble element stiffness vector
            vector_slice<c_vector<double, 2*ELEMENT_DIM+2> > rBElem_slice_V  (rBElem, slice (0, 2, num_elem_nodes));
            vector_slice<c_vector<double, 2*ELEMENT_DIM+2> > rBElem_slice_Phi(rBElem, slice (1, 2, num_elem_nodes));
            
            rBElem_slice_V   += wJ*( (Am*Cm*Vm/mDt - Am*I_ionic - I_intra_stim) * basis_func );
            rBElem_slice_Phi += wJ*( -I_extra_stim * basis_func );
        }
    } 
    
    
    
    void AssembleSystem(Vec currentSolution, double currentTime)
    {
        // Allow the PDE to set up anything necessary for the assembly of the
        // solution (ie solve the ODEs)
        mpBidomainPde->PrepareForAssembleSystem(currentSolution, currentTime);
                
        if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->ZeroRhsVector();
        } 
        else 
        {
            mpAssembledLinearSystem->ZeroLinearSystem();
        }
     
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::MeshIterator iter =
            mpMesh->GetElementIteratorBegin();
 

        // assumes elements all have same number of nodes
        const int num_elem_nodes = iter->GetNumNodes();
        

        c_matrix<double, 2*ELEMENT_DIM+2, 2*ELEMENT_DIM+2> a_elem;
        c_vector<double, 2*ELEMENT_DIM+2> b_elem;
        
        while (iter != mpMesh->GetElementIteratorEnd())
        {
            const Element<ELEMENT_DIM, SPACE_DIM> &element = *iter;

            AssembleOnElement(element, a_elem, b_elem);
         
            for (int i=0; i<num_elem_nodes; i++)
            {
                int node1 = element.GetNodeGlobalIndex(i);
                
                if (!mMatrixIsAssembled)
                {
                    for (int j=0; j<num_elem_nodes; j++)
                    {
                        int node2 = element.GetNodeGlobalIndex(j);
 
                        mpAssembledLinearSystem->AddToMatrixElement(2*node1,   2*node2,   a_elem(2*i,   2*j));
                        mpAssembledLinearSystem->AddToMatrixElement(2*node1+1, 2*node2,   a_elem(2*i+1, 2*j));
                        mpAssembledLinearSystem->AddToMatrixElement(2*node1,   2*node2+1, a_elem(2*i,   2*j+1));
                        mpAssembledLinearSystem->AddToMatrixElement(2*node1+1, 2*node2+1, a_elem(2*i+1, 2*j+1));
                    }
                }
                
                mpAssembledLinearSystem->AddToRhsVectorElement(2*node1,   b_elem(2*i));
                mpAssembledLinearSystem->AddToRhsVectorElement(2*node1+1, b_elem(2*i+1));
            }
            iter++;
        }

        
        if (!mMatrixIsAssembled)
        {
            //create null space for matrix and pass to linear system
            Vec nullbasis[1];  
            unsigned lo, hi;
            mpBidomainPde->GetOwnershipRange(lo, hi);
            VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*mpMesh->GetNumNodes(), &nullbasis[0]);  
            double* p_nullbasis;
            VecGetArray(nullbasis[0], &p_nullbasis); 

            for (unsigned global_index=lo; global_index<hi; global_index++)
            {
                unsigned local_index = global_index - lo;
                p_nullbasis[2*local_index  ] = 0;
                p_nullbasis[2*local_index+1] = 1;
            }
            VecRestoreArray(nullbasis[0], &p_nullbasis);      
            VecAssemblyBegin(nullbasis[0]);
            VecAssemblyEnd(nullbasis[0]);
            
            mpAssembledLinearSystem->SetNullBasis(nullbasis, 1);
            
            
            VecDestroy(nullbasis[0]);
        }
    



        if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->AssembleRhsVector();
        } 
        else 
        {
            mpAssembledLinearSystem->AssembleIntermediateLinearSystem();
        }
        
        // no need to loop over surface elems to apply neumann bcs here, because the 
        // neumann bcs are zero
                
        // no need to apply dirichlet boundary conditions for phi either (?-check with Jon) 
        // (code left here for debugging/comparison) - 
        /*
        //       for(int node_num = 468; node_num<507; node_num++)       
        int node_num = 0;
        {
            if(!mMatrixIsAssembled)
            {
                mpAssembledLinearSystem->ZeroMatrixRow   ( 2*node_num + 1 );
                mpAssembledLinearSystem->SetMatrixElement( 2*node_num + 1, 2*node_num + 1,1);
            }
            mpAssembledLinearSystem->SetRhsVectorElement( 2*node_num + 1, 0);
        }
        */
        
        if (mMatrixIsAssembled)
        {
            mpAssembledLinearSystem->AssembleRhsVector();
        } 
        else 
        {
            mpAssembledLinearSystem->AssembleFinalLinearSystem();
        }
        
        mMatrixIsAssembled = true;


        // write matrix and rhs vector
        //mpAssembledLinearSystem->WriteLinearSystem("mat.txt", "rhs.txt");       
    }
   
    
public:
    BidomainDg0Assembler(BidomainPde<SPACE_DIM>* pBidomainPde, 
                         ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh, 
                         AbstractLinearSolver* pLinearSolver,
                         int numQuadPoints = 2)
      : AbstractAssembler<ELEMENT_DIM,SPACE_DIM>(numQuadPoints),
	mpAssembledLinearSystem(NULL)
    {
        mpBidomainPde = pBidomainPde;
        mpMesh = pMesh;
        mpSolver = pLinearSolver;

        mMatrixIsAssembled = false;
    }
    
    
    ~BidomainDg0Assembler()
    {
        if (mpAssembledLinearSystem != NULL)
        {
            delete mpAssembledLinearSystem;
            mpAssembledLinearSystem = NULL;
        }
    }
   
    
    // identical to SetTimes in SimpleDg0ParabolicAssembler
    void SetTimes(double Tstart, double Tend, double dT)
    {
        mTstart = Tstart;
        mTend   = Tend;
        mDt     = dT;
        
        assert(mTstart < mTend);
        assert(mDt > 0);
        assert(mDt <= mTend - mTstart + 1e-12);
    
        mTimesSet = true;
    }
    
    
    // identical to SetInitialCondition in SimpleDg0ParabolicAssembler
    void SetInitialCondition(Vec initCondition)
    {
        /// \todo: check initCondition is the correct size, & do the same in other assemblers
        mInitialCondition = initCondition;
        mInitialConditionSet = true;

        if (mpAssembledLinearSystem == NULL)
        {
            mpAssembledLinearSystem = new LinearSystem(initCondition);
        }
    }
    

    // very similar to Solve in SimpleDg0ParabolicAssembler.    
    // (this doesn't take in any arguments however, and also the AssembleSystem
    // method in this assembler class does not call Solve() on the linear system,
    // instead it is done here).
    Vec Solve()
    {
        assert(mTimesSet);
        assert(mInitialConditionSet);
        
        double t = mTstart;
 
        Vec current_solution = mInitialCondition;
        Vec next_solution;
        while( t < mTend - 1e-10 )
        {
            AssembleSystem(current_solution, t);

            // std::cout << "solving at t="<< t << "\n" <<std::flush;
            next_solution = mpAssembledLinearSystem->Solve(mpSolver);
                    
            t += mDt;

            if (current_solution != mInitialCondition)
            {
                VecDestroy(current_solution);
            }
            current_solution = next_solution;
        }   
        
        return current_solution;
    }      
};




#endif /*_BIDOMAINDG0ASSEMBLER_HPP_*/
