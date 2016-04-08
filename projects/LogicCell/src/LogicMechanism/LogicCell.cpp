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


#include "CellActions.hpp"
#include "LogicCell.hpp"

#include "WildTypeCellMutationState.hpp"
#include "DummyCellCycleModel.hpp"


LogicCell::LogicCell(boost::shared_ptr<AbstractCellProperty> pDummyMut, bool archiving, CellPropertyCollection prop):
        Cell( pDummyMut, new DummyCellCycleModel(), NULL, archiving, prop)
{
};


LogicCell::~LogicCell(){
   for(int i = 0; i < (int)cellLogicModules.size(); i++){
      if(cellLogicModules[i]!=NULL) {
         delete cellLogicModules[i];
      }
   }
}


void LogicCell::SetCellCycleModel(AbstractCellCycleModel* pCellCycleModel){
   EXCEPTION("SetCellCycleModel is not permitted; LogicCells always have a DummyCellCycleModel");
};

void LogicCell::SetMutationState(boost::shared_ptr<AbstractCellProperty> pMutationState){
   EXCEPTION("SetMutationState is not permitted; LogicCells always have a DummyMutationState");
};

void LogicCell::SetBirthTime(double birthTime){
   mBirthTime = birthTime;
};

double LogicCell::GetBirthTime() const{
   return mBirthTime;
};


void LogicCell::InitialiseCellCycleModel(){};


double LogicCell::GetAge() const{
   return SimulationTime::Instance()->GetTime() - mBirthTime;
}


bool LogicCell::ReadyToDivide(){

   assert(!IsDead());
   if (mUndergoingApoptosis)
   {
      return false;
   }

   UpdateLogic();
   return mCanDivide;
};



CellPtr LogicCell::Divide(){

   assert(!IsDead());
   assert(mCanDivide);
   mCanDivide = false;

   // Create daughter cell
   LogicCell* p_new_cell = MakeNewCell();
   p_new_cell->SetApoptosisTime(mApoptosisTime);
   p_new_cell->mBirthTime = SimulationTime::Instance()->GetTime();


   return boost::shared_ptr<Cell>(p_new_cell);
};



LogicCell* LogicCell::MakeNewCell(){

   CellPropertyCollection daughterProp = rGetCellPropertyCollection();
   daughterProp.RemoveProperty<CellId>();
   daughterProp.RemoveProperty<CellData>();

   LogicCell* daughterCell = new LogicCell(GetMutationState(), false, daughterProp);
   daughterCell->cellLogicModules.resize(cellLogicModules.size(), NULL);

   for(int i = 0; i < (int) cellLogicModules.size(); i++){
      if(cellLogicModules[i] != NULL) {
         AbstractCellLogic *daughterLogic = cellLogicModules[i]->Divide(daughterCell);
         daughterLogic->SetCategory(cellLogicModules[i]->GetCategory());
         daughterLogic->DumpState();
         daughterCell->SetLogicById(i,daughterLogic);
      }
   }

   return daughterCell;
}



void LogicCell::UpdateLogic() {

   ClearPendingActions();

   for (int i = 0; i < (int) cellLogicModules.size(); i++) {

      if (cellLogicModules[i] != NULL) {
         cellLogicModules[i]->Update();
      }
   }

   ProcessPendingActions();
}



void LogicCell::ProcessPendingActions(){

   for(int i=0; i < (int)pendingOutboundActions.size(); i++){

      int action = pendingOutboundActions[i];
      if(action == CellActions::CallForDivision){
         this->mCanDivide = true;
      }
   }
}


std::vector<int> LogicCell::GetPendingActions(){
   return pendingOutboundActions;
};


void LogicCell::ClearPendingActions(){
   pendingOutboundActions.clear();
};


void LogicCell::SetPendingAction(int actionCode){
   pendingOutboundActions.push_back(actionCode);
}


void LogicCell::SetLogicById(int positionInVector, AbstractCellLogic* newlogic){
   cellLogicModules[positionInVector] = newlogic;
};
