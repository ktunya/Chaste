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

#include "PCLDUFactorisation.hpp"
#include "Exception.hpp"

PCLDUFactorisation::PCLDUFactorisation(KSP& rKspObject)
{    
    PCLDUFactorisationCreate(rKspObject);    
    PCLDUFactorisationSetUp();
}
    
PCLDUFactorisation::~PCLDUFactorisation()
{
    MatDestroy(mPCContext.A11_matrix_subblock);
    MatDestroy(mPCContext.A22_matrix_subblock);
    MatDestroy(mPCContext.B_matrix_subblock);
    
    PCDestroy(mPCContext.PC_amg_A11);
    PCDestroy(mPCContext.PC_amg_A22);
    
    VecDestroy(mPCContext.x11);
    VecDestroy(mPCContext.y11);

    VecDestroy(mPCContext.x22);
    VecDestroy(mPCContext.y22);

    VecDestroy(mPCContext.z);
    VecDestroy(mPCContext.temp);
}

void PCLDUFactorisation::PCLDUFactorisationCreate(KSP& rKspObject)
{
    KSPGetPC(rKspObject, &mPetscPCObject);        
    
    Mat system_matrix, dummy;
    MatStructure flag;    
    KSPGetOperators(rKspObject, &system_matrix, &dummy, &flag); 
    
    PetscInt num_rows, num_columns;
    MatGetSize(system_matrix, &num_rows, &num_columns);
    assert(num_rows==num_columns);  
    
    PCSetType(mPetscPCObject, PCSHELL);
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PCShellSetApply(mPetscPCObject, PCLDUFactorisationApply, (void*) &mPCContext);
#else
    // Register PC context so it gets passed to PCBlockDiagonalApply
    PCShellSetContext(mPetscPCObject, &mPCContext);
    // Register call-back function
    PCShellSetApply(mPetscPCObject, PCLDUFactorisationApply);
#endif
        
    // Get matrix sublock A11
    IS A11_rows, A11_columns;    
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 2, &A11_rows);
    ISCreateStride(PETSC_COMM_WORLD, num_columns/2, 0, 2, &A11_columns);
    
    MatGetSubMatrix(system_matrix, A11_rows, A11_columns, PETSC_DECIDE, MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);

    ISDestroy(A11_rows);
    ISDestroy(A11_columns);    

    // Get matrix sublock A22
    IS A22_rows, A22_columns;    
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 1, 2, &A22_rows);
    ISCreateStride(PETSC_COMM_WORLD, num_columns/2, 1, 2, &A22_columns);
    
    MatGetSubMatrix(system_matrix, A22_rows, A22_columns, PETSC_DECIDE, MAT_INITIAL_MATRIX, &mPCContext.A22_matrix_subblock);

    ISDestroy(A22_rows);
    ISDestroy(A22_columns);    

    // Get matrix sublock B (the upper triangular one)
    IS B_rows, B_columns;    
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 2, &B_rows);
    ISCreateStride(PETSC_COMM_WORLD, num_columns/2, 1, 2, &B_columns);
    
    MatGetSubMatrix(system_matrix, B_rows, B_columns, PETSC_DECIDE, MAT_INITIAL_MATRIX, &mPCContext.B_matrix_subblock);

    ISDestroy(B_rows);
    ISDestroy(B_columns);
    
    // Allocate memory
    mPCContext.x11 = PetscTools::CreateVec(num_rows/2);
    mPCContext.x22 = PetscTools::CreateVec(num_rows/2);
    mPCContext.y11 = PetscTools::CreateVec(num_rows/2);
    mPCContext.y22 = PetscTools::CreateVec(num_rows/2);       
    mPCContext.z = PetscTools::CreateVec(num_rows/2);
    mPCContext.temp = PetscTools::CreateVec(num_rows/2);
}

void PCLDUFactorisation::PCLDUFactorisationSetUp()
{
    // These options will get read by PCSetFromOptions
//    PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");        
//    PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");            
//    PetscOptionsSetValue("-pc_hypre_type", "boomeramg");        

    // Set up amg preconditioner for block A11
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A11));
    PCSetOperators(mPCContext.PC_amg_A11, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);

////////
//    PCSetType(mPCContext.PC_amg_A11, PCHYPRE);
////////
    PCSetType(mPCContext.PC_amg_A11, PCKSP);
    KSP ksp1;
    PCKSPGetKSP(mPCContext.PC_amg_A11,&ksp1);
    KSPSetType(ksp1, KSPCG);
    KSPSetTolerances(ksp1, 0.1, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);        

    PC prec1;
    KSPGetPC(ksp1, &prec1);
    PCSetType(prec1, PCBJACOBI);
    PCSetFromOptions(prec1);
    PCSetOperators(prec1, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
    PCSetUp(prec1);

    KSPSetFromOptions(ksp1);
    KSPSetUp(ksp1);
////////

    PCSetFromOptions(mPCContext.PC_amg_A11);
    PCSetUp(mPCContext.PC_amg_A11);


    // Set up amg preconditioner for block A22
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A22));
    PCSetOperators(mPCContext.PC_amg_A22, mPCContext.A22_matrix_subblock, mPCContext.A22_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);

////////
//    PCSetType(mPCContext.PC_amg_A22, PCHYPRE);
////////
    PCSetType(mPCContext.PC_amg_A22, PCKSP);
    KSP ksp2;
    PCKSPGetKSP(mPCContext.PC_amg_A22,&ksp2);
    KSPSetType(ksp2, KSPCG);
    KSPSetTolerances(ksp2, 0.1, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);    

    PC prec2;
    KSPGetPC(ksp2, &prec2);
    PCSetType(prec2, PCBJACOBI);
    PCSetFromOptions(prec2);
    PCSetOperators(prec2, mPCContext.A22_matrix_subblock, mPCContext.A22_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);    
    PCSetUp(prec2);

    KSPSetFromOptions(ksp2);
    KSPSetUp(ksp2);
////////

    PCSetFromOptions(mPCContext.PC_amg_A22);
    PCSetUp(mPCContext.PC_amg_A22);        
}

PetscErrorCode PCLDUFactorisationApply(void* pc_context, Vec x, Vec y)
{
    /// \todo refactoring: create a method for scattering and another for reversing
    
    // Cast the pointer to a PC context to our defined type
    PCLDUFactorisation::PCLDUFactorisationContext* block_diag_context = (PCLDUFactorisation::PCLDUFactorisationContext*) pc_context;
    assert(block_diag_context!=NULL); 
    
    /////////////////////
    PetscInt num_rows;
    VecGetSize(x, &num_rows);

    IS A11_rows;    
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 2, &A11_rows);    
    
    VecScatter A11_scatter_ctx;
    VecScatterCreate(x, A11_rows, block_diag_context->x11, PETSC_NULL, &A11_scatter_ctx);

//PETSc-3.x.x or PETSc-2.3.3 
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(A11_scatter_ctx, x, block_diag_context->x11, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(A11_scatter_ctx, x, block_diag_context->x11, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x11, INSERT_VALUES, SCATTER_FORWARD, A11_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x11, INSERT_VALUES, SCATTER_FORWARD, A11_scatter_ctx);
#endif    

    IS A22_rows;
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 1, 2, &A22_rows);

    VecScatter A22_scatter_ctx;
    VecScatterCreate(x, A22_rows, block_diag_context->x22, PETSC_NULL, &A22_scatter_ctx);

//PETSc-3.x.x or PETSc-2.3.3 
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(A22_scatter_ctx, x, block_diag_context->x22, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(A22_scatter_ctx, x, block_diag_context->x22, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x22, INSERT_VALUES, SCATTER_FORWARD, A22_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x22, INSERT_VALUES, SCATTER_FORWARD, A22_scatter_ctx);
#endif    

    ////////////////////
    
    //z  = inv(A11)*x1
    PCApply(block_diag_context->PC_amg_A11, block_diag_context->x11, block_diag_context->z);
    double minus_one = -1.0;
    //y2 = inv(A22)*(x2 - B*z)
    MatMult(block_diag_context->B_matrix_subblock,block_diag_context->z,block_diag_context->temp); //temp = B*z    
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    VecAYPX(&minus_one, block_diag_context->x22, block_diag_context->temp); // temp <-- x2 - temp
#else
    VecAYPX(block_diag_context->temp, minus_one, block_diag_context->x22); // temp <-- x2 - temp
#endif
    PCApply(block_diag_context->PC_amg_A22, block_diag_context->temp, block_diag_context->y22); // y2 = inv(A22)*temp
    
    //y1 = z - inv(A11)(B*y2)
    MatMult(block_diag_context->B_matrix_subblock,block_diag_context->y22,block_diag_context->temp); //temp = B*y2 
    ///\todo Are these lines in the correct order Miguel?
    PCApply(block_diag_context->PC_amg_A11, block_diag_context->temp, block_diag_context->y11); // y1 = inv(A11)*temp
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    VecAYPX(&minus_one, block_diag_context->z, block_diag_context->y11); // y1 <-- z - y1

#else
    VecAYPX(block_diag_context->y11, minus_one, block_diag_context->z); // y1 <-- z - y1
#endif
                
    ////////////////////

//PETSc-3.x.x or PETSc-2.3.3 
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(A11_scatter_ctx, block_diag_context->y11, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(A11_scatter_ctx, block_diag_context->y11, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y11, y, INSERT_VALUES, SCATTER_REVERSE, A11_scatter_ctx);
    VecScatterEnd(block_diag_context->y11, y, INSERT_VALUES, SCATTER_REVERSE, A11_scatter_ctx);
#endif    

//PETSc-3.x.x or PETSc-2.3.3 
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(A22_scatter_ctx, block_diag_context->y22, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(A22_scatter_ctx, block_diag_context->y22, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y22, y, INSERT_VALUES, SCATTER_REVERSE, A22_scatter_ctx);
    VecScatterEnd(block_diag_context->y22, y, INSERT_VALUES, SCATTER_REVERSE, A22_scatter_ctx);
#endif    
    
    ////////////////////

    ISDestroy(A11_rows);
    ISDestroy(A22_rows);
        
    VecScatterDestroy(A11_scatter_ctx);
    VecScatterDestroy(A22_scatter_ctx);    
    
    return 0;
}    
