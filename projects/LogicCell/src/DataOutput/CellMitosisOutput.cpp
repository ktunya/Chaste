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

#include "CellMitosisOutput.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "LogicTypes.hpp"

//Constructor, initialises sampling interval and sets output file to null
template<unsigned DIM>
CellMitosisOutput<DIM>::CellMitosisOutput(int interval, double startRow, double endRow, double spacing)
    : AbstractCellBasedSimulationModifier<DIM>(),
      OutputFile(NULL),
      mInterval(interval),
      mStartRow(startRow),
      mEndRow(endRow),
      mSpacing(spacing)   
{
}


//Empty destructor 
template<unsigned DIM>
CellMitosisOutput<DIM>::~CellMitosisOutput(){}


//Getter methods for private members
template<unsigned DIM>
int CellMitosisOutput<DIM>::GetInterval() const
{
  return mInterval;
};
template<unsigned DIM>
double CellMitosisOutput<DIM>::GetStartX() const{
  return mStartRow;
};
template<unsigned DIM>
double CellMitosisOutput<DIM>::GetEndX() const{
  return mEndRow;
};
template<unsigned DIM>
double CellMitosisOutput<DIM>::GetSpacing() const{
  return mSpacing;
};



//Open an output file MitosisData.txt in the simulation directory
template<unsigned DIM>
void CellMitosisOutput<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
  OutputFileHandler rOutputFileHandler(outputDirectory, false);
  OutputFile = rOutputFileHandler.OpenOutputFile("MitosisData.txt");
}


//At each timestep, if the time is a sampling time, loop through all cells and compile data.
//Output that data to file.
template<unsigned DIM>
void CellMitosisOutput<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

  

  //If it's a sampling time, start gathering data
  if(SimulationTime::Instance()->GetTimeStepsElapsed() % GetInterval() == 0){

    int nRows = (int)((mEndRow - mStartRow)/mSpacing)+2;

    std::vector<int> cellCount = std::vector<int>(nRows,0);
    std::vector<int> mitosisCount = std::vector<int>(nRows,0);

    for (typename AbstractCellPopulation<DIM,DIM>::RealCellsIterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End(); ++cell_iter)
    { 

      Node<DIM>* node = rCellPopulation.GetNode(rCellPopulation.GetLocationIndexUsingCell(*cell_iter));
      c_vector<double, DIM> loc = node->rGetLocation();

      int rowIndex = (int)((loc[0]-mStartRow)/mSpacing);
      cellCount[rowIndex] = cellCount[rowIndex]+1;      

      int state = cell_iter->GetCellData()->GetItem("CellCycle");
      if(state == CellCycle::M){
        mitosisCount[rowIndex] = mitosisCount[rowIndex]+1;
      }

    }

    //Write data
    *OutputFile << SimulationTime::Instance()->GetTime() << "\t";
    for(int i=0; i<cellCount.size(); i++){
      *OutputFile << cellCount[i] << "\t";
    }
    for(int i=0; i<mitosisCount.size(); i++){
      *OutputFile << mitosisCount[i];
      if(i!=(mitosisCount.size()-1)){
        *OutputFile << "\t"; 
      } 
    }
    *OutputFile << "\n";

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
void CellMitosisOutput<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<SampleMeiDataEveryXTimesteps>" << GetInterval() << "</SampleMeiDataEveryXTimesteps>\n";
  // Call method on direct parent class
  AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile); 
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellMitosisOutput<1>;
template class CellMitosisOutput<2>;
template class CellMitosisOutput<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellMitosisOutput)
