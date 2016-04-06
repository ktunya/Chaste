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

#include "DTCSignal.hpp"
#include "EnvironmentalSignalTypes.hpp"
#include "LogicCell.hpp"

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
DTCSignal<ELEMENT_DIM, SPACE_DIM>::DTCSignal():
  AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>(),
  linearDecayBegins(DOUBLE_UNSET),
  signalIsZero(DOUBLE_UNSET)
{
};


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void DTCSignal<ELEMENT_DIM, SPACE_DIM>::SetSignalProperties(double decayBegins, double signalZero){
  linearDecayBegins = decayBegins;
  signalIsZero = signalZero;
};


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void DTCSignal<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    for(typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::RealCellsIterator it = rCellPopulation.Begin();
        it!=rCellPopulation.End(); ++it){

        c_vector<double, SPACE_DIM> position = rCellPopulation.GetLocationOfCellCentre(*it);
        double signalStrength;
        if(position[0]<linearDecayBegins){
          signalStrength = 1;
        }else if(position[0]>=linearDecayBegins && position[0]<signalIsZero){
          signalStrength = 1- (position[0]-linearDecayBegins)/(signalIsZero - linearDecayBegins);
        }else{
          signalStrength = 0;
        }

        it->GetCellData()->SetItem("DTCSignal", signalStrength);
        dynamic_cast<LogicCell*>((*it).get())->template sendEnvironmentalSignal<MyDTCSignal>(signalStrength);
    }
}


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void DTCSignal<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    for(typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::RealCellsIterator it = rCellPopulation.Begin();
        it!=rCellPopulation.End(); ++it){

        c_vector<double, SPACE_DIM> position = rCellPopulation.GetLocationOfCellCentre(*it);
        double signalStrength;
        if(position[0]<linearDecayBegins){
          signalStrength = 1;
        }else if(position[0]>=linearDecayBegins && position[0]<signalIsZero){
          signalStrength = 1- (position[0]-linearDecayBegins)/(signalIsZero - linearDecayBegins);
        }else{
          signalStrength = 0;
        }

        it->GetCellData()->SetItem("DTCSignal", signalStrength);
        dynamic_cast<LogicCell*>((*it).get())->template sendEnvironmentalSignal<MyDTCSignal>(signalStrength);

    }
}


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void DTCSignal<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile){
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
};

// Explicit instantiation
template class DTCSignal<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DTCSignal)