/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "RK4NumericalMethod.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
RK4NumericalMethod<ELEMENT_DIM,SPACE_DIM> :: RK4NumericalMethod():
  AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>()
{	
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
RK4NumericalMethod<ELEMENT_DIM,SPACE_DIM>::~RK4NumericalMethod(){
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void RK4NumericalMethod<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt){

    if(this->mNonEulerSteppersEnabled){

      std::vector<c_vector<double, SPACE_DIM> > K1 = this->ComputeAndSaveForcesInclDamping();

      // Update node positions by +(dt K1/2) and compute K2
      int index = 0;
      for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->pCellPopulation->rGetMesh().GetNodeIteratorBegin();
       node_iter != this->pCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
      {
        c_vector<double, SPACE_DIM> newLocation = node_iter->rGetLocation() + dt * K1[index]/2.0;
        this->SafeNodePositionUpdate(node_iter->GetIndex(), newLocation);
      }
      std::vector< c_vector<double, SPACE_DIM> > K2 = this->ComputeAndSaveForcesInclDamping(); 


      // Update node positions by +(dt K2-K1 /2) and compute K3
      index = 0;
      for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->pCellPopulation->rGetMesh().GetNodeIteratorBegin();
           node_iter != this->pCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
      {
        c_vector<double, SPACE_DIM> newLocation = node_iter->rGetLocation() + dt * (K2[index] - K1[index])/2.0; //revert, then update
        this->SafeNodePositionUpdate(node_iter->GetIndex(), newLocation);
      }              
      std::vector< c_vector<double, SPACE_DIM> > K3 = this->ComputeAndSaveForcesInclDamping(); 


      // Update node positions by +(dt K3-K2 /2) and compute K4
      index = 0;
      for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->pCellPopulation->rGetMesh().GetNodeIteratorBegin();
           node_iter != this->pCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
      {
        double damping = this->pCellPopulation->GetDampingConstant(node_iter->GetIndex());
        c_vector<double, SPACE_DIM> newLocation = node_iter->rGetLocation() + dt * (K3[index] - K2[index]/2.0); //revert, then update
        this->SafeNodePositionUpdate(node_iter->GetIndex(), newLocation);
      }
      std::vector< c_vector<double, SPACE_DIM> > K4 = this->ComputeAndSaveForcesInclDamping(); 


      // Final position update, according to the RK4 numerical method
      index = 0;
      for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->pCellPopulation->rGetMesh().GetNodeIteratorBegin();
       node_iter != this->pCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
      {
        c_vector<double, SPACE_DIM> effectiveForce = (K1[index] + 2*K2[index] + 2*K3[index] + K4[index])/6.0;
        c_vector<double, SPACE_DIM> oldLocation = node_iter->rGetLocation() - dt * K3[index]; //revert
        c_vector<double, SPACE_DIM> finalDisplacement =  dt * effectiveForce;

        this->DetectStepSizeExceptions(node_iter->GetIndex(), finalDisplacement, dt);
        
        c_vector<double, SPACE_DIM> newLocation = oldLocation + finalDisplacement;     
        this->SafeNodePositionUpdate(node_iter->GetIndex(), newLocation);

        // Ensure that each nodes holds an accurate applied force value, 
        // incase it's accessed by some other class
        double damping = this->pCellPopulation->GetDampingConstant(node_iter->GetIndex());
        node_iter->ClearAppliedForce();
        c_vector<double, SPACE_DIM> force = effectiveForce*damping;
        node_iter->AddAppliedForceContribution(force);
      }

    }else{

      // If this type of cell population does not support the new numerical methods, delegate 
      // updating node positions to the population itself.
      // Basically only applies to NodeBasedCellPopulationWithBuskeUpdates.     
      this->ComputeAndSaveForcesInclDamping();
      this->pCellPopulation->UpdateNodeLocations(dt);    

    }
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void RK4NumericalMethod<ELEMENT_DIM, SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile){
  // Nothing yet
};


///////// Explicit instantiation
template class RK4NumericalMethod<1,1>;
template class RK4NumericalMethod<1,2>;
template class RK4NumericalMethod<2,2>;
template class RK4NumericalMethod<1,3>;
template class RK4NumericalMethod<2,3>;
template class RK4NumericalMethod<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RK4NumericalMethod)