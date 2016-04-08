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

#ifndef ABSTRACTCELLLOGIC_H
#define ABSTRACTCELLLOGIC_H

#include <string>

class LogicCell;

class AbstractCellLogic{

protected:

    // Current state of this logic. Often an enum for readability. 
    // E.g. the state of a CellCycle logic could be CellCycle::G1 (i.e. 0)
    int state;

    // A descriptive name for this type of logic module. Will be used whenever
    // data is output for visualisation. So for example, in Paraview there'll be
    // a property called "category" associated with each cell containing the state info.
    std::string category;

    // Pointer to the cell that contains this logic module.
    LogicCell* owningCell;

public:

    // Constructor. Takes a pointer to the owning cell, and an int indicating the 
    // starting state.
    AbstractCellLogic(LogicCell* inputCell, int initialState);
    virtual ~AbstractCellLogic();

    // Pure virtual method. Tells the logic how to update each timestep.
    virtual void Update() = 0;

    // Produces a daughter logic module, with its starting state determined by that of the
    // parent.
    virtual AbstractCellLogic* Divide(LogicCell* owningCell) = 0;

    // Returns the current state of the logic module
    int GetState();

    // Sets the descriptive name for this type of logic module (std::string category)
    void SetCategory(std::string inputString);

    // Gets the descriptive name for this type of logic module (std::string category)
    std::string GetCategory();

    // Outputs data for visualisation.
    virtual void DumpState();

};

#endif //ABSTRACTCELLLOGIC_H
