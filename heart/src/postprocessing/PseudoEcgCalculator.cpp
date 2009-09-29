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

#include "PseudoEcgCalculator.hpp"
#include "PetscTools.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> ::GetIntegrand(ChastePoint<SPACE_DIM>& rX,
                                c_vector<double,PROBLEM_DIM>& rU,
                                c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU)
{
    double denominator = 0;
    for (unsigned i = 0; i < SPACE_DIM; i++)
    {
        denominator = denominator + (rX[i] - mrX[i])*(rX[i] - mrX[i]); 
    }
    
    c_vector<double,SPACE_DIM> grad_one_over_r;   
    for (unsigned j = 0; j < SPACE_DIM; j++)
    {
        grad_one_over_r[j] = - (rX[j] - mrX[j])*pow( (1/denominator) , 1.5); 
    }
    
    double integrand = 0;
    for (unsigned k = 0; k < SPACE_DIM; k++)
    {
        integrand = integrand + rGradU(k, k) * grad_one_over_r[k]; 
    }
    
    return -mDiffusionCoefficient*integrand;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> ::PseudoEcgCalculator (TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh, 
                                                                                 ChastePoint<SPACE_DIM>& rX,
                                                                                 std::string directory, 
                                                                                 std::string hdf5File, 
                                                                                 bool makeAbsolute)
                                      : mrMesh(rMesh),
                                        mrX(rX)
                                      
{
    mpDataReader = new Hdf5DataReader(directory, hdf5File, makeAbsolute);
    mNumberOfNodes = mpDataReader->GetNumberOfRows();
    mNumTimeSteps = mpDataReader->GetUnlimitedDimensionValues().size();
    mDiffusionCoefficient = 1.0;
    //check that the hdf file was generated by simulations from the same mesh
    assert(mNumberOfNodes == mrMesh.GetNumNodes());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~PseudoEcgCalculator()
{
    delete mpDataReader;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetDiffusionCoefficient(double diffusionCoefficient)
{
    assert(diffusionCoefficient>=0);
    mDiffusionCoefficient = diffusionCoefficient;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
std::vector<double> PseudoEcgCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputePseudoEcg ()
{
    Vec SolutionAtOneTimestep = PetscTools::CreateVec(mNumberOfNodes);
    std::vector<double> pseudo_ecg;
    for (unsigned i = 0; i < mNumTimeSteps; i++)
    {
        mpDataReader->GetVariableOverNodes(SolutionAtOneTimestep, "V" , i);
        double pseudo_ecg_at_one_timestep = Calculate(mrMesh, SolutionAtOneTimestep);
        pseudo_ecg.push_back(pseudo_ecg_at_one_timestep);
    }
    VecDestroy(SolutionAtOneTimestep);
    return pseudo_ecg;
}
/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class PseudoEcgCalculator<1,1,1>;
template class PseudoEcgCalculator<1,2,1>;
template class PseudoEcgCalculator<1,3,1>;
template class PseudoEcgCalculator<1,2,2>;
template class PseudoEcgCalculator<2,3,1>;
template class PseudoEcgCalculator<3,3,1>;
