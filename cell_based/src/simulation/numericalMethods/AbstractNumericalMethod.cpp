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

#include "AbstractNumericalMethod.hpp"
#include "AbstractMesh.hpp"
#include "StepSizeException.hpp"
#include "Warnings.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "CellBasedEventHandler.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::AbstractNumericalMethod()
:pCellPopulation(NULL),
 pForceCollection(NULL)
{	
    /* pCellPopulation and pForceCollection are initialized by OffLatticeSimulation in
    its constructor. For now we set some safe defaults for the other member variables */

    mNonEulerSteppersEnabled = false;
    mGhostNodeForcesEnabled = true;
    mParticleForcesEnabled = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::~AbstractNumericalMethod(){	
};



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetCellPopulation(AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>* pPopulation){
    
    if(!pCellPopulation){
        pCellPopulation = pPopulation;
    }else{
        EXCEPTION("The cell population referred to by a numerical method should not be reset");
    }    

    // Set other member variables according to the type of the cell population
    if(dynamic_cast<NodeBasedCellPopulationWithBuskeUpdate<SPACE_DIM>*>(pCellPopulation)){
        mNonEulerSteppersEnabled = false;
        WARNING("Non-Euler steppers are not yet implemented for NodeBasedCellPopulationWithBuskeUpdate");
    }else{
        mNonEulerSteppersEnabled = true;
    }

    if(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>*>(pCellPopulation)){
        mGhostNodeForcesEnabled = true;
    }else{
        mGhostNodeForcesEnabled = false;
    }

    if(dynamic_cast<NodeBasedCellPopulationWithParticles<SPACE_DIM>*>(pCellPopulation)){
        mParticleForcesEnabled = true;
    }else{
        mParticleForcesEnabled = false;
    }
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SetForceCollection(std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >* pForces){
    
    if(!pForceCollection){
        pForceCollection = pForces;
    }else{
        EXCEPTION("The force collection referred to by a numerical method should not be reset");
    }    
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::ComputeAndSaveForces(){
	
	CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);

    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }
    
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = pForceCollection->begin();
        iter != pForceCollection->end(); ++iter)
    {
        (*iter)->AddForceContribution(*pCellPopulation);
    }

    // Here deal with forces on non-cell nodes (ghosts and particles)
    if(mGhostNodeForcesEnabled){
        //dynamic_cast<MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>*>(&rCellPopulation)->ApplyGhostForces();
    }
    if(mParticleForcesEnabled){
        //dynamic_cast<NodeBasedCellPopulationWithParticles<SPACE_DIM>*>(&rCellPopulation)->ApplyParticleForces();
    }
    WARN_ONCE_ONLY("REENABLE SPECIAL CASE NODE FORCES (AbstractNumericalMethod)");

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

    // Store the forces in a vector
	std::vector<c_vector<double, SPACE_DIM> > forcesAsVector;
	forcesAsVector.reserve(pCellPopulation->GetNumNodes());

    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter)
    {
        double damping = pCellPopulation->GetDampingConstant(node_iter->GetIndex());
        forcesAsVector.push_back(node_iter->rGetAppliedForce()/damping); 
    }

    return forcesAsVector;
};



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::SaveCurrentLocations(){

	std::vector<c_vector<double, SPACE_DIM> > currentLocations;
	currentLocations.reserve(pCellPopulation->GetNumNodes());

	int index = 0;
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd(); 
         ++index, ++node_iter)
    {
		currentLocations[index] = node_iter->rGetLocation();	
	}

	return currentLocations;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::HandleStepSizeExceptions(c_vector<double,SPACE_DIM>* displacement, double dt, unsigned nodeIndex){
    
    try{
        WARN_ONCE_ONLY("REENABLE STEP SIZE EXCEPTIONS (AbstractNumericalMethod)");
        //pCellPopulation->FindAndAddressStepSizeExceptions(displacement, dt, nodeIndex);
    
    }catch(StepSizeException* e){

        if(!(e->isTerminal)){
            WARN_ONCE_ONLY(e->what());
        }else{
            throw e;
        }
    }   
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodInfo(out_stream& rParamsFile){

    std::string numerical_method_type = GetIdentifier();

    *rParamsFile << "\t\t<" << numerical_method_type << ">\n";
    OutputNumericalMethodParameters(rParamsFile);
    *rParamsFile << "\t\t</" << numerical_method_type << ">\n";
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile){
    // May be implemented in child classes
};


///////// Explicit instantiation
template class AbstractNumericalMethod<1,1>;
template class AbstractNumericalMethod<1,2>;
template class AbstractNumericalMethod<2,2>;
template class AbstractNumericalMethod<1,3>;
template class AbstractNumericalMethod<2,3>;
template class AbstractNumericalMethod<3,3>;