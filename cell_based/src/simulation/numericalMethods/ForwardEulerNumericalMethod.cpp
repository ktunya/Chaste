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

#include "ForwardEulerNumericalMethod.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
ForwardEulerNumericalMethod<ELEMENT_DIM,SPACE_DIM>::ForwardEulerNumericalMethod()
:AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>()
{	
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
ForwardEulerNumericalMethod<ELEMENT_DIM,SPACE_DIM>::~ForwardEulerNumericalMethod(){
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void ForwardEulerNumericalMethod<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt){

  if(this->mNonEulerSteppersEnabled){

    // Apply forces to each cell, and save a vector of net forces F
    std::vector<c_vector<double, SPACE_DIM> > F = this->ComputeAndSaveForcesInclDamping();

    int index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->pCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->pCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      // Get the current node location and calculate the new location according to forward Euler 
      c_vector<double, SPACE_DIM> oldLocation = node_iter->rGetLocation();
      c_vector<double, SPACE_DIM> displacement = dt * F[index];

      // In the vertex-based case, displacement may be scaled if the cellrearrangement threshold is exceeded 
      this->DetectStepSizeExceptions(node_iter->GetIndex(), displacement, dt);  

      c_vector<double, SPACE_DIM> newLocation = oldLocation + displacement;
      this->SafeNodePositionUpdate(node_iter->GetIndex(), newLocation);        
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
void ForwardEulerNumericalMethod<ELEMENT_DIM, SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile){
  // Nothing here yet
};


///////// Explicit instantiation
template class ForwardEulerNumericalMethod<1,1>;
template class ForwardEulerNumericalMethod<1,2>;
template class ForwardEulerNumericalMethod<2,2>;
template class ForwardEulerNumericalMethod<1,3>;
template class ForwardEulerNumericalMethod<2,3>;
template class ForwardEulerNumericalMethod<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ForwardEulerNumericalMethod)