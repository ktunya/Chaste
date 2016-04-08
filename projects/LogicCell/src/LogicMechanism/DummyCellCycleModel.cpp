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

#include "DummyCellCycleModel.hpp"
#include <stdio.h>

DummyCellCycleModel::DummyCellCycleModel():
AbstractCellCycleModel(){

};


DummyCellCycleModel::~DummyCellCycleModel(){};


void DummyCellCycleModel::Initialise(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. Initialise is ignored." << std::endl;
};


void DummyCellCycleModel::InitialiseDaughterCell(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. InitialiseDaughterCell is ignored." << std::endl;
};


void DummyCellCycleModel::SetBirthTime(double birthTime){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. SetBirthTime is ignored." << std::endl;
};


double DummyCellCycleModel::GetBirthTime() const{
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetBirthTime is ignored." << std::endl;
   return 0.0;
};


double DummyCellCycleModel::GetAge(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetAge is ignored." << std::endl;
   return 0.0;
};


bool DummyCellCycleModel::ReadyToDivide(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. ReadyToDivide is ignored." << std::endl;
   return false;
};


void DummyCellCycleModel::UpdateCellCyclePhase(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. UpdateCellCyclePhase is ignored." << std::endl;
};


AbstractCellCycleModel* DummyCellCycleModel::CreateCellCycleModel(){

   DummyCellCycleModel* p_model = new DummyCellCycleModel();
   return p_model;
};


void DummyCellCycleModel::ResetForDivision(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. ResetForDivision is ignored." << std::endl;
};


CellCyclePhase DummyCellCycleModel::GetCurrentCellCyclePhase(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetCurrentCellCyclePhase is ignored." << std::endl;
   return G_ZERO_PHASE;
};


double DummyCellCycleModel::GetG1Duration(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetG1Duration is ignored." << std::endl;
   return 0.0;
};


double DummyCellCycleModel::GetStemCellG1Duration(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetStemCellG1Duration is ignored." << std::endl;
   return 0.0;
};


double DummyCellCycleModel::GetTransitCellG1Duration(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetTransitCellG1Duration is ignored." << std::endl;
   return 0.0;
};


double DummyCellCycleModel::GetSG2MDuration(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetSG2MDuration is ignored." << std::endl;
   return 0.0;
};


double DummyCellCycleModel::GetSDuration(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetSDuration is ignored." << std::endl;
   return 0.0;
};


double DummyCellCycleModel::GetG2Duration(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetG2Duration is ignored." << std::endl;
   return 0.0;
};


double DummyCellCycleModel::GetMDuration(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetMDuration is ignored." << std::endl;
   return 0.0;
};


void DummyCellCycleModel::SetStemCellG1Duration(double stemCellG1Duration){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. SetStemCellG1Duration is ignored." << std::endl;
};


void DummyCellCycleModel::SetTransitCellG1Duration(double transitCellG1Duration){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. SetTransitCellG1Duration is ignored." << std::endl;
};


void DummyCellCycleModel::SetSDuration(double sDuration){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. SetSDuration is ignored." << std::endl;
};


void DummyCellCycleModel::SetG2Duration(double g2Duration){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. SetG2Duration is ignored." << std::endl;
};


void DummyCellCycleModel::SetMDuration(double mDuration){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. SetMDuration is ignored." << std::endl;
};


double DummyCellCycleModel::GetAverageTransitCellCycleTime(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetAverageTransitCellCycleTime is ignored." << std::endl;
   return 0.0;
};


double DummyCellCycleModel::GetAverageStemCellCycleTime(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetAverageStemCellCycleTime is ignored." << std::endl;
   return 0.0;
};


bool DummyCellCycleModel::CanCellTerminallyDifferentiate(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. CanCellTerminallyDifferentiate is ignored." << std::endl;
   return false;
};


double DummyCellCycleModel::GetMinimumGapDuration(){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. GetMinimumGapDuration is ignored." << std::endl;
   return 0.0;
};


void DummyCellCycleModel::SetMinimumGapDuration(double minimumGapDuration){
   std::cout << "WARNING: Logic cells use a CellCycleLogic not a CellCycleModel. SetMinimumGapDuration is ignored." << std::endl;
};


void DummyCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile){

   *rParamsFile << "\t\t\t<CCM_WARNING>" << " DUMMY CELL CYCLE MODEL IS UNUSED; PLEASE SEE CELLCYCLE LOGIC MODULE " << "</CCM_WARNING>\n";
   AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DummyCellCycleModel)