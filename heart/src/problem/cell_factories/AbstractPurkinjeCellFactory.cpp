/*

Copyright (c) 2005-2012, University of Oxford.
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


#include "AbstractPurkinjeCellFactory.hpp"
#include "PurkinjeVentricularJunctionStimulus.hpp"
#include "MultiStimulus.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::AbstractPurkinjeCellFactory()
    : AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>(),
      mpMixedDimensionMesh(NULL)
{
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::CreateJunction(AbstractCardiacCell* pPurkinjeCell,
                                                                        AbstractCardiacCell* pCardiacCell,
                                                                        double resistance)
{
    // Create the junction stimuli, and associate them with the cells
    boost::shared_ptr<PurkinjeVentricularJunctionStimulus> p_pvj_ventricular_stim(new PurkinjeVentricularJunctionStimulus(resistance));
    boost::shared_ptr<PurkinjeVentricularJunctionStimulus> p_pvj_purkinje_stim(new PurkinjeVentricularJunctionStimulus(resistance));
    p_pvj_purkinje_stim->SetAppliedToPurkinjeCellModel();
    p_pvj_ventricular_stim->SetVentricularCellModel(pCardiacCell);
    p_pvj_ventricular_stim->SetPurkinjeCellModel(pPurkinjeCell);
    p_pvj_purkinje_stim->SetVentricularCellModel(pCardiacCell);
    p_pvj_purkinje_stim->SetPurkinjeCellModel(pPurkinjeCell);

    // Create new combined stimuli which add the junction stimuli to those already in the cells
    boost::shared_ptr<MultiStimulus> p_multi_stim_ventricular(new MultiStimulus);
    p_multi_stim_ventricular->AddStimulus(p_pvj_ventricular_stim);
    p_multi_stim_ventricular->AddStimulus(pCardiacCell->GetStimulusFunction());
    pCardiacCell->SetStimulusFunction(p_multi_stim_ventricular);

    boost::shared_ptr<MultiStimulus> p_multi_stim_purkinje(new MultiStimulus);
    p_multi_stim_purkinje->AddStimulus(p_pvj_purkinje_stim);
    p_multi_stim_purkinje->AddStimulus(pPurkinjeCell->GetStimulusFunction());
    pPurkinjeCell->SetStimulusFunction(p_multi_stim_purkinje);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    mpMixedDimensionMesh = dynamic_cast<MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>*>(pMesh);
    if (mpMixedDimensionMesh ==NULL)
    {
        EXCEPTION("AbstractPurkinjeCellFactory must take a MixedDimensionMesh");
    }
    mLocalPurkinjeNodes.clear();
    for (typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::CableElementIterator iter = mpMixedDimensionMesh->GetCableElementIteratorBegin();
          iter != mpMixedDimensionMesh->GetCableElementIteratorEnd();
          ++iter)
    {
        mLocalPurkinjeNodes.insert((*iter)->GetNodeGlobalIndex(0u));
        mLocalPurkinjeNodes.insert((*iter)->GetNodeGlobalIndex(1u));
    }
    AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(pMesh);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCell*  AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::CreatePurkinjeCellForNode(
        unsigned nodeIndex,
        AbstractCardiacCell* pCardiacCell)
{
    if (mLocalPurkinjeNodes.count(nodeIndex)>0)
    {
        return CreatePurkinjeCellForTissueNode(nodeIndex, pCardiacCell);
    }
    else
    {
        return new FakeBathCell(this->mpSolver, this->mpZeroStimulus);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::GetMixedDimensionMesh()
{
    if (mpMixedDimensionMesh == NULL)
    {
        EXCEPTION("The mixed dimension mesh object has not been set in the cell factory");
    }
    return mpMixedDimensionMesh;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractPurkinjeCellFactory<1,1>;
template class AbstractPurkinjeCellFactory<2,2>;
template class AbstractPurkinjeCellFactory<3,3>;
template class AbstractPurkinjeCellFactory<1,2>;
template class AbstractPurkinjeCellFactory<1,3>;
