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
#include "SimpleAssymmetricStemDivision.hpp"
#include "Exception.hpp"
#include "LogicTypes.hpp"
#include "LogicCell.hpp"


SimpleAssymmetricStemDivision::SimpleAssymmetricStemDivision(LogicCell* inputCell, int initialState, int initialTADivisions):
        AbstractCellLogic(inputCell, initialState),
        nTADivisions(initialTADivisions)
{

}


void SimpleAssymmetricStemDivision::update(){

   if(state == Differentiation::D){
      return;
   }


   if(nTADivisions > 2 && owningCell->getState<CellCycle>() == CellCycle::G1){

      state = Differentiation::D;

   }else{

      double nicheSignal = owningCell->getEnvironmentalSignal<NicheSignal>();
      //owningCell->GetCellData()->SetItem("NicheSignal", nicheSignal);
      std::cout << "NicheSignal: " << nicheSignal << std::endl;

      if(nicheSignal > 5){
         state = Differentiation::TA;
      }

   }

   dumpState();
}


AbstractCellLogic* SimpleAssymmetricStemDivision::divide(LogicCell* daughterCell){


   if(state == Differentiation::S){

      AbstractCellLogic* daughterLogic =  new SimpleAssymmetricStemDivision(daughterCell, Differentiation::TA, 0);
      dumpState();
      return daughterLogic;

   }else if(state == Differentiation::TA) {

      AbstractCellLogic* daughterLogic =  new SimpleAssymmetricStemDivision(daughterCell, Differentiation::TA, nTADivisions+1);
      nTADivisions = nTADivisions+1;
      dumpState();
      return daughterLogic;

   }else if(state == Differentiation::D){
      EXCEPTION("Differentiated cells cannot divide");
   }

   EXCEPTION("Invalid state in SimpleAssymmetricStemDivision");
   return NULL;
}


void SimpleAssymmetricStemDivision::dumpState(){
   owningCell->GetCellData()->SetItem(category, state);
}


