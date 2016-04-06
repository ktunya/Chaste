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

#include "CellGrowthWithDiffDependence.hpp"

#include "LogicTypes.hpp"
#include "LogicCell.hpp"


CellGrowthWithDiffDependence::CellGrowthWithDiffDependence(LogicCell* inputCell, 
                                                           int initialState,
                                                           double initialRadius):
  AbstractCellLogic(inputCell, initialState),
  cellRadius(initialRadius)
{
   radialGrowthPerHour = DOUBLE_UNSET;
   maximumRadiusUndiff = DOUBLE_UNSET;
   maximumRadiusDiff = DOUBLE_UNSET;
}


void CellGrowthWithDiffDependence::setHourlyRadialGrowthRate( double rate ){
   radialGrowthPerHour = rate;
};
void CellGrowthWithDiffDependence::setMaxRadiusUndiffd( double maxRadUndiffd ){
   maximumRadiusUndiff = maxRadUndiffd;
};
void CellGrowthWithDiffDependence::setMaxRadiusDiffd( double maxRadDiffd ){
   maximumRadiusDiff = maxRadDiffd;
};


void CellGrowthWithDiffDependence::update(){

    if(radialGrowthPerHour == DOUBLE_UNSET){
      EXCEPTION("Please call setHourlyRadialGrowthRate on CellGrowthWithDiffDependence");
    }
    if(maximumRadiusDiff == DOUBLE_UNSET){
      EXCEPTION("Please call setMaxRadiusDiffd on CellGrowthWithDiffDependence");
    }
    if(maximumRadiusUndiff == DOUBLE_UNSET){
      EXCEPTION("Please call setMaxRadiusUndiffd on CellGrowthWithDiffDependence");
    }

    double maxRadius;
    if(owningCell->getState<CellCycle>() == CellCycle::G0){
      maxRadius = maximumRadiusDiff;
    }else{
      maxRadius = maximumRadiusUndiff;
    }

    if(cellRadius >= maxRadius) {
        state = Growth::STATIC;
    }else{
        state = Growth::GROWING;
    }

    if(state == Growth::GROWING) {
        cellRadius += radialGrowthPerHour * SimulationTime::Instance()->GetTimeStep();
    }  

    dumpState();
}


AbstractCellLogic* CellGrowthWithDiffDependence::divide(LogicCell* daughterCell){

    CellGrowthWithDiffDependence* daughterLogic =  new CellGrowthWithDiffDependence(daughterCell, Growth::GROWING, cellRadius/pow(2,0.33333333));
    
    daughterLogic->setHourlyRadialGrowthRate(radialGrowthPerHour);
    daughterLogic->setMaxRadiusDiffd(maximumRadiusDiff);
    daughterLogic->setMaxRadiusUndiffd(maximumRadiusUndiff);

    cellRadius = cellRadius/pow(2,0.33333333);
    dumpState();
    daughterLogic->dumpState();

    return daughterLogic;
}


void CellGrowthWithDiffDependence::dumpState(){
    owningCell->GetCellData()->SetItem(category, state);
    owningCell->GetCellData()->SetItem("Radius", cellRadius);
}