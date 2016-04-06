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

#ifndef LOGICTYPES_H
#define LOGICTYPES_H

#include <iostream>
using namespace std;

class AllBaseLogicTypes
{
protected:
    static int numberOfLogicTypesInUse;  //initialized to zero
};

template <typename T>
class BaseLogicType : public AllBaseLogicTypes
{
protected:
    static int logicTypeId;

public:
    static int id()
    {
        if (logicTypeId == -1)
        {
            logicTypeId = AllBaseLogicTypes::numberOfLogicTypesInUse++;
        }
        return logicTypeId;
    }
};

template <typename T>
int BaseLogicType<T>::logicTypeId = -1; //Initialise ID of each logic type to -1


//A macro for declaring new types of cell logic. Shields users from having to
//think about the curiously recursive template pattern:
#define CLASS_LOGIC_TYPE(x) class x : public BaseLogicType<x>


CLASS_LOGIC_TYPE(CellCycle){
public:

    enum CellCyclePhase {G1, S, G2, M, MeiG1, MeiS, G0};

};


CLASS_LOGIC_TYPE(Differentiation){

public:

    enum DifferentionState {S, TA, D};

};


CLASS_LOGIC_TYPE(Growth){

public:

    enum GrowthState {STATIC, GROWING};

};


#endif //LOGICTYPES_H