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

#include "ZoneOutputOutput.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "LogicTypes.hpp"

//Constructor, sets output file to null
template<unsigned DIM>
ZoneOutputOutput<DIM>::ZoneOutputOutput(int stepsPerHour)
    : AbstractCellBasedSimulationModifier<DIM>(),
      OutputFile(NULL),
      StepsPerHour(stepsPerHour)
{}


//Empty destructor 
template<unsigned DIM>
ZoneOutputOutput<DIM>::~ZoneOutputOutput(){}


//Open an output file OutputData.txt in the simulation directory
template<unsigned DIM>
void ZoneOutputOutput<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
  OutputFileHandler rOutputFileHandler(outputDirectory, false);
  OutputFile = rOutputFileHandler.OpenOutputFile("OutputData.txt");
}


//Getter methods for private members
template<unsigned DIM>
int ZoneOutputOutput<DIM>::GetStepsPerHour() const
{
  return StepsPerHour;
};



//Compile data and output hourly
template<unsigned DIM>
void ZoneOutputOutput<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

  int stepsElapsed = SimulationTime::Instance()->GetTimeStepsElapsed();

  if( stepsElapsed % StepsPerHour == 0){

    int NewMeioticCellsThisHour = 0;

    for (typename AbstractCellPopulation<DIM,DIM>::RealCellsIterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End(); ++cell_iter)
    { 

      int state = cell_iter->GetCellData()->GetItem("CellCycle");
    
      if(state == CellCycle::G0){

        double TimeAtHim3Labelling = cell_iter->GetCellData()->GetItem("TimeAtHim3Labelling");

        if(TimeAtHim3Labelling > (SimulationTime::Instance()->GetTime() - 1) ){

          NewMeioticCellsThisHour++;
        }

      }

    }

    //Write data
    *OutputFile << SimulationTime::Instance()->GetTime() << "\t" 
                << NewMeioticCellsThisHour << "\n";

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
void ZoneOutputOutput<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  // Call method on direct parent class
  AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile); 
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ZoneOutputOutput<1>;
template class ZoneOutputOutput<2>;
template class ZoneOutputOutput<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ZoneOutputOutput)
