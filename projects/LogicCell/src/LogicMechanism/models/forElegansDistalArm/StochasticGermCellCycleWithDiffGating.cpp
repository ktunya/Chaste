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

#include "StochasticGermCellCycleWithDiffGating.hpp"
#include "CellActions.hpp"
#include "RandomNumberGenerator.hpp"

StochasticGermCellCycleWithDiffGating::StochasticGermCellCycleWithDiffGating(LogicCell* inputCell, 
                                                                             int initialState, 
                                                                             double initialTimeInPhase_Hrs):
  AbstractCellLogic(inputCell, initialState),
  timeInPhase(initialTimeInPhase_Hrs)
{
    phaseDurations = std::vector<double>();
    stochasticity = DOUBLE_UNSET;
    currentPhaseDuration = DOUBLE_UNSET;
    gate = -1;
}


void StochasticGermCellCycleWithDiffGating::setPhaseDurations( std::vector<double> durations ){
    phaseDurations = durations;
};
void StochasticGermCellCycleWithDiffGating::setStochasticity( double stochasticityParam ){
    stochasticity = stochasticityParam;
};
void StochasticGermCellCycleWithDiffGating::setGate( int gatingState ){
    gate = gatingState;
};


void StochasticGermCellCycleWithDiffGating::update(){

  if(gate == -1){
    EXCEPTION("Please call setGatingState on StochasticGermCellCycleWithDiffGating and provide an integer indicating from which phase entry into meiosis is allowed to occur; G1 = 0, S = 1, G2 = 2, M = 3");
  }
  if(currentPhaseDuration == DOUBLE_UNSET){
    getNewPhaseDuration(state);
  }


  //Increment time spent in this phase
  timeInPhase += SimulationTime::Instance()->GetTimeStep();


  //Check whether a differentiation signal has been sent, and the gating state reached. 
  //If so enter the meiosis pathway
  if( state == gate && owningCell->getState<Differentiation>() == Differentiation::D){        
    state = CellCycle::MeiG1;
    timeInPhase = 0.0;
    getNewPhaseDuration(state);
  } 


  //If the current cell cycle phase is complete...
  if (timeInPhase >= currentPhaseDuration) {

    //For mitotic cells, cycle the state and reset timeInPhase / currentPhaseDuration
    if(state < 4){   
      state++;
      state = state%4;
      timeInPhase = 0.0;
      getNewPhaseDuration(state);

    }

    //For meiotic pathway cells, send G1->MeioticS and MeioticS->G0
      if(state == CellCycle::MeiG1){
         state = CellCycle::MeiS;
         timeInPhase = 0.0;
         getNewPhaseDuration(state); 
      }else if(state == CellCycle::MeiS){
         owningCell->GetCellData()->SetItem("TimeAtHim3Labelling", SimulationTime::Instance()->GetTime());
         state = CellCycle::G0;
         currentPhaseDuration = std::numeric_limits<double>::max();
      }
  }

  //If M phase was just entered, request a cell division from the mechanics simulation 
  if (state == CellCycle::M && timeInPhase == 0.0) {
     owningCell->SetPendingAction(CellActions::CallForDivision);
  }

  //output state for visualisation
  dumpState();
}


void StochasticGermCellCycleWithDiffGating::getNewPhaseDuration(int st){

    if(phaseDurations.size()!=6){
    EXCEPTION("Please call setPhaseDurations on StochasticGermCellCycleWithDiffGating and provide lengths in hours for G1, S, G2, M, \"Meiotic G1\" and Meiotic S.");
    }
    if(stochasticity == DOUBLE_UNSET){
    EXCEPTION("Please call setStochasticity on StochasticGermCellCycleWithDiffGating and provide a number between 0 and 1 indicating the desired standard deviation as a proportion of the mean");
    }

    double meanDuration = phaseDurations[st];

    if(st == CellCycle::G1 || st == CellCycle::G2){
      double sd = stochasticity * meanDuration;
      RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
      currentPhaseDuration = p_gen->NormalRandomDeviate(meanDuration, sd);
    }else{
      currentPhaseDuration = meanDuration;
    }
}


AbstractCellLogic* StochasticGermCellCycleWithDiffGating::divide(LogicCell* daughterCell){

  StochasticGermCellCycleWithDiffGating* daughterLogic =  new StochasticGermCellCycleWithDiffGating(daughterCell, state, timeInPhase);
  
  daughterLogic->setPhaseDurations(phaseDurations);
  daughterLogic->setStochasticity(stochasticity);
  daughterLogic->setGate(gate);

  //Make sure correct data is output for visualisation
  dumpState();
  daughterLogic->dumpState(); 
  
  return daughterLogic;
}


void StochasticGermCellCycleWithDiffGating::dumpState(){
  owningCell->GetCellData()->SetItem(category, state);
}