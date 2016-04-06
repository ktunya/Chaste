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

#include "EnvironmentalSignalTypes.hpp"
#include "PositionSignal.hpp"
#include "LogicCell.hpp"

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
PositionSignal<ELEMENT_DIM, SPACE_DIM>::PositionSignal():
  AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>(),
  axis(-1),
  zero(DOUBLE_UNSET)
{
};


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void PositionSignal<ELEMENT_DIM, SPACE_DIM>::SetAxisOfInterest(int axisChoice){
  axis = axisChoice;
};

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void PositionSignal<ELEMENT_DIM, SPACE_DIM>::SetMeasurementStartPoint(double z){
  zero = z;
};


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void PositionSignal<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{

  if(zero == DOUBLE_UNSET){
    std::cout << "Please call SetMeasurementStartPoint on PositionSignal and specify the location to measure positions from" << std::endl;
    EXCEPTION("PositionSignal exception");
  }

  for(typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::RealCellsIterator it = rCellPopulation.Begin();
      it!=rCellPopulation.End(); ++it){

      c_vector<double, SPACE_DIM> position = rCellPopulation.GetLocationOfCellCentre(*it);
      double signalStrength;
      if(axis==0){
        signalStrength = position[0]-zero;
      }else if(axis==1){
        signalStrength = position[1]-zero;
      }else if(axis==2){
        signalStrength = position[2]-zero;
      }else{
        std::cout << "Please call SetAxisOfInterest on PositionSignal and specify an axis between 0 and 2" << std::endl;
        EXCEPTION("PositionSignal exception");
      }

      it->GetCellData()->SetItem("PositionSignal", signalStrength);
      dynamic_cast<LogicCell*>((*it).get())->template sendEnvironmentalSignal<MyPositionSignal>(signalStrength);
  }
}


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void PositionSignal<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{

  if(zero == DOUBLE_UNSET){
    std::cout << "Please call SetMeasurementStartPoint on PositionSignal and specify the location to measure positions from" << std::endl;
    EXCEPTION("PositionSignal exception");
  }

  for(typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::RealCellsIterator it = rCellPopulation.Begin();
      it!=rCellPopulation.End(); ++it){

      c_vector<double, SPACE_DIM> position = rCellPopulation.GetLocationOfCellCentre(*it);
      double signalStrength;
      if(axis==0){
        signalStrength = position[0]-zero;
      }else if(axis==1){
        signalStrength = position[1]-zero;
      }else if(axis==2){
        signalStrength = position[2]-zero;
      }else{
        std::cout << "Please call SetAxisOfInterest on PositionSignal and specify an axis between 0 and 2" << std::endl;
        EXCEPTION("PositionSignal exception");
      }

      it->GetCellData()->SetItem("PositionSignal", signalStrength);
      dynamic_cast<LogicCell*>((*it).get())->template sendEnvironmentalSignal<MyPositionSignal>(signalStrength);
  }
}


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void PositionSignal<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile){
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
};

// Explicit instantiation
template class PositionSignal<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PositionSignal)
