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

#include "FixedDivisionTimingsCellCycleModel.hpp"
#include "Warnings.hpp"


FixedDivisionTimingsCellCycleModel::FixedDivisionTimingsCellCycleModel(double timeBetweenDivisions):
targetTimeBetweenDivisions(timeBetweenDivisions)
{
	WARN_ONCE_ONLY("WARNING: A FixedDivisionTimingsCellCycleModel manually divides cells every T hours,\
	allowing variability in division timings to be controlled when comparing different simulations.\
	It does not update or use cell cycle phases, and cell cycle phase getters and setters will\
	return irrelevant or unset values.");
}	

FixedDivisionTimingsCellCycleModel::~FixedDivisionTimingsCellCycleModel()
{
}


void FixedDivisionTimingsCellCycleModel::UpdateCellCyclePhase(){

	this->mReadyToDivide = false;
	double currentTime = SimulationTime::Instance()->GetTime();	

	if(currentTime != 0){
	
		int cycleNumber = (int)(currentTime / targetTimeBetweenDivisions);
		double currentRemainder = fmin( fabs(currentTime - (cycleNumber*targetTimeBetweenDivisions)),  
			                            fabs(((cycleNumber+1)*targetTimeBetweenDivisions) - currentTime) );
	
		double threshold = SimulationTime::Instance()->GetTimeStep()/2.0;
	
		if(currentRemainder == threshold && (currentTime - (cycleNumber*targetTimeBetweenDivisions)) < 0 ){
			this->mReadyToDivide = true;
		}
	
		if(currentRemainder < threshold){
			this->mReadyToDivide = true;
		}

	}
};


AbstractCellCycleModel* FixedDivisionTimingsCellCycleModel::CreateCellCycleModel(){
	return new FixedDivisionTimingsCellCycleModel(targetTimeBetweenDivisions);
};


void FixedDivisionTimingsCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile){
	*rParamsFile << "\t\t\t<IntervalBetweenDivisions>" << targetTimeBetweenDivisions << "</IntervalBetweenDivisions>\n";
};


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedDivisionTimingsCellCycleModel)