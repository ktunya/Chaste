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

#include "GenerationalDifferentiationModel.hpp"
#include "Exception.hpp"
#include "LogicTypes.hpp"
#include "EnvironmentalSignalTypes.hpp"
#include "LogicCell.hpp"

GenerationalDifferentiationModel::GenerationalDifferentiationModel(LogicCell* inputCell, int initialState):
  AbstractCellLogic(inputCell, initialState){ 

    maxDivisions = -1;
    nDivisions = 0;
};


void GenerationalDifferentiationModel::setMaxDivisions( int divisions ){
   maxDivisions = divisions;
};

void GenerationalDifferentiationModel::setNDivisions( int divisions ){
   nDivisions = divisions;
};


void GenerationalDifferentiationModel::update(){

   if(maxDivisions == -1){
      EXCEPTION("Please call setnDivisions on GenerationalDifferentiationModel");
   }

   dumpState();
}



AbstractCellLogic* GenerationalDifferentiationModel::divide(LogicCell* daughterCell){

   GenerationalDifferentiationModel* daughterLogic;
   if(state == Differentiation::S){
      //Assymmetric division 
      daughterLogic = new GenerationalDifferentiationModel(daughterCell, Differentiation::TA);
      daughterLogic -> setMaxDivisions(maxDivisions);
      daughterLogic -> setNDivisions(0);

   }else if(state==Differentiation::TA){

      if(nDivisions < maxDivisions){
        //TA cell that has not exhausted its proliferative capacity. Becomes 2 mew TAs.
        daughterLogic = new GenerationalDifferentiationModel(daughterCell, Differentiation::TA);
        daughterLogic -> setMaxDivisions(maxDivisions);
        daughterLogic -> setNDivisions(nDivisions+1);
        nDivisions++;
      }else{
        //TA cell that HAS exhausted its proliferative capacity
        daughterLogic = new GenerationalDifferentiationModel(daughterCell, Differentiation::D);
        daughterLogic -> setMaxDivisions(maxDivisions);
        daughterLogic -> setNDivisions(nDivisions);
        state = Differentiation::D;
      }

   }else if(state == Differentiation::D){
        //Extra division by a germ cell that is preparing to differentiate, but is awaiting
        //G1 before it can do so.
        daughterLogic = new GenerationalDifferentiationModel(daughterCell, Differentiation::D);
        daughterLogic -> setMaxDivisions(maxDivisions);
        daughterLogic -> setNDivisions(nDivisions);
   }

   return daughterLogic;

};


void GenerationalDifferentiationModel::dumpState(){
    owningCell->GetCellData()->SetItem(category, state);
}

