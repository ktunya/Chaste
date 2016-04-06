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

#include "CellCyclePhaseOutput.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "LogicTypes.hpp"

//Constructor, initialises sampling interval and sets output file to null
template<unsigned DIM>
CellCyclePhaseOutput<DIM>::CellCyclePhaseOutput(int interval)
    : AbstractCellBasedSimulationModifier<DIM>(),
      OutputFile(NULL),
      mInterval(interval)
{}


//Empty destructor 
template<unsigned DIM>
CellCyclePhaseOutput<DIM>::~CellCyclePhaseOutput(){}


//Getter methods for private members
template<unsigned DIM>
int CellCyclePhaseOutput<DIM>::GetInterval() const
{
  return mInterval;
};


//Open an output file PhaseData.txt in the simulation directory
template<unsigned DIM>
void CellCyclePhaseOutput<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
  OutputFileHandler rOutputFileHandler(outputDirectory, false);
  OutputFile = rOutputFileHandler.OpenOutputFile("PhaseData.txt");
}


//At each timestep, if the time is a sampling time, loop through all cells and compile data.
//Output that data to file.
template<unsigned DIM>
void CellCyclePhaseOutput<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

  //If it's a sampling time, start gathering data
  if(SimulationTime::Instance()->GetTimeStepsElapsed() % GetInterval() == 0){

    int G1Count = 0;
    int SCount = 0;
    int G2Count = 0;
    int MCount = 0;
    int MeiG1Count = 0;
    int MeiSCount = 0;
    int DiffCount = 0;

    double firstMeiotic = 250;
    double lastMitotic = 0;

    for (typename AbstractCellPopulation<DIM,DIM>::RealCellsIterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End(); ++cell_iter)
    { 

      Node<DIM>* node = rCellPopulation.GetNode(rCellPopulation.GetLocationIndexUsingCell(*cell_iter));
      c_vector<double, DIM> loc = node->rGetLocation();
      double xPos = loc[0];

      int state = cell_iter->GetCellData()->GetItem("CellCycle");
      if(state == CellCycle::G1){
        if(xPos > lastMitotic){
          lastMitotic = xPos;
        } 
        G1Count++;
      }else if(state == CellCycle::S){
        if(xPos > lastMitotic){
          lastMitotic = xPos;
        } 
        SCount++;
      }else if(state == CellCycle::G2){
        if(xPos > lastMitotic){
          lastMitotic = xPos;
        } 
        G2Count++;
      }else if(state == CellCycle::M){
        if(xPos > lastMitotic){
          lastMitotic = xPos;
        } 
        MCount++;
      }else if(state == CellCycle::MeiG1){
        if(xPos > lastMitotic){
          lastMitotic = xPos;
        } 
        MeiG1Count++;  
      }else if(state == CellCycle::MeiS){
        if(xPos > lastMitotic){
          lastMitotic = xPos;
        } 
        MeiSCount++;
      }else if(state == CellCycle::G0){
        if(xPos < firstMeiotic){
          firstMeiotic = xPos;
        } 
        DiffCount++;
      }else{
        std::cout << "INVALID CELL CYCLE STATE " << state << " DETECTED\n";
      }

    }

    if( G1Count+SCount+G2Count+MCount+MeiG1Count+MeiSCount > 750 ){
      std::cout << "Max prolif cells exceeded" << std::endl;
      EXCEPTION("Max prolif cells exceeded");
    }

    //Write data
    *OutputFile << SimulationTime::Instance()->GetTime() << "\t" 
                << G1Count << "\t" 
                << SCount << "\t" 
                << G2Count << "\t" 
                << MCount << "\t"
                << MeiG1Count << "\t"
                << MeiSCount << "\t"
                << DiffCount << "\t"
                << lastMitotic << "\t"
                << firstMeiotic << "\n";

    //Flush the output file to record data as soon as possible
    OutputFile->flush();
  }

  //If the simulation is finished, close the output file.
  if(SimulationTime::Instance()->IsFinished()){
    OutputFile->close();
  }

}


//Output this class's parameters to a log file
template<unsigned DIM>
void CellCyclePhaseOutput<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<SamplePhaseDataEveryXTimesteps>" << GetInterval() << "</SamplePhaseDataEveryXTimesteps>\n";
  // Call method on direct parent class
  AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile); 
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellCyclePhaseOutput<1>;
template class CellCyclePhaseOutput<2>;
template class CellCyclePhaseOutput<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellCyclePhaseOutput)
