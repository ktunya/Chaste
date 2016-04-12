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

#include "CellCycleTestWithDivisionMechanism.hpp"
#include "CellActions.hpp"
#include "Warnings.hpp"

CellCycleTestWithDivisionMechanism::CellCycleTestWithDivisionMechanism(LogicCell* inputCell, int initialState, double startingTime, double phaseDuration):
  AbstractCellLogic(inputCell, initialState),
  phaseStartTime(startingTime),
  phaseDuration(phaseDuration)
{
}


void CellCycleTestWithDivisionMechanism::Update(){
    
    int oldState = state;

    double timeInPhase;
    if(SimulationTime::Instance()->IsStartTimeSetUp()){
      timeInPhase = SimulationTime::Instance()->GetTime()-phaseStartTime;
    }else{
      // We're running inside a test where cells are updated manually rather than
      // via a simulation object. This allows the cell cycle to work without a 
      // SimulationTime being available.
      phaseStartTime--;
      timeInPhase = -phaseStartTime;
    }

    if (timeInPhase >= phaseDuration) {
       if(SimulationTime::Instance()->IsStartTimeSetUp()){
          phaseStartTime = SimulationTime::Instance()->GetTime();
       }else{
          phaseStartTime = 0.0;
       }
       state++;
       state = state % 4;
    }

    if (state == CellCycle::M && oldState == CellCycle::G2) {
      owningCell->SetPendingAction(CellActions::CallForDivision);
    }

    DumpState();
}

AbstractCellLogic* CellCycleTestWithDivisionMechanism::Divide(LogicCell* daughterCell){

    AbstractCellLogic* daughterLogic;
    if(SimulationTime::Instance()->IsStartTimeSetUp()){
      daughterLogic =  new CellCycleTestWithDivisionMechanism(daughterCell, 
                                                              state, 
                                                              SimulationTime::Instance()->GetTime(),
                                                              phaseDuration);
    }else{
      daughterLogic =  new CellCycleTestWithDivisionMechanism(daughterCell, 
                                                              state, 
                                                              0,
                                                              phaseDuration);
    }

    DumpState();
    daughterLogic->DumpState();

    return daughterLogic;
}
