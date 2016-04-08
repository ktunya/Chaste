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

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"
#include "Warnings.hpp"

#include "LogicCell.hpp"
#include "AbstractCellLogic.hpp"
#include "LogicTypes.hpp"
#include "EnvironmentalSignalTypes.hpp"
#include "CellActions.hpp"
#include "DummyMutationState.hpp"

#include "FixedDurationTestCellCycle.hpp"
#include "SimpleAssymmetricStemDivision.hpp"
#include "DifferentiateIfSignalAboveThresh.hpp"
#include "CellCycleTestWithDivisionMechanism.hpp"



#ifndef TESTLOGICCELL_HPP
#define TESTLOGICCELL_HPP

class TestLogicCell : public AbstractCellBasedTestSuite {


public:

    void TestInitialisation(){

        MAKE_PTR(DummyMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        FixedDurationTestCellCycle* cycleLogic = new FixedDurationTestCellCycle(newCell, CellCycle::G1, 0.0, 2.0);
        newCell->SetLogic<CellCycle>(cycleLogic);

        int state = newCell->GetState<CellCycle>();
        TS_ASSERT_EQUALS(state, CellCycle::G1);

        for(int i=0; i<5; i++){
            newCell->UpdateLogic();
        }

        state = newCell->GetState<CellCycle>();
        TS_ASSERT_EQUALS(state, CellCycle::G2);

        delete newCell;
    }

    
    void TestSymmetricDivision(){

        MAKE_PTR(DummyMutationState, pDummyMut);
        LogicCell* parentCell = new LogicCell(pDummyMut);

        FixedDurationTestCellCycle* cycleLogic = new FixedDurationTestCellCycle(parentCell, CellCycle::G1, 0.0, 2.0);
        parentCell->SetLogic<CellCycle>(cycleLogic);

        LogicCell* daughterCell = parentCell->MakeNewCell();

        TS_ASSERT_EQUALS(parentCell->GetState<CellCycle>(), daughterCell->GetState<CellCycle>());

        delete parentCell;
        delete daughterCell; 
    }


    void TestAssymmetricDivision(){

        MAKE_PTR(DummyMutationState, pDummyMut);
        LogicCell* parentCell = new LogicCell(pDummyMut);

        FixedDurationTestCellCycle* cycleLogic = new FixedDurationTestCellCycle(parentCell, CellCycle::M, 0.0, 2.0);
        parentCell->SetLogic<CellCycle>(cycleLogic);
        SimpleAssymmetricStemDivision* diffLogic = new SimpleAssymmetricStemDivision(parentCell, Differentiation::S, 0);
        parentCell->SetLogic<Differentiation>(diffLogic);

        LogicCell* daughterCell = parentCell->MakeNewCell();

        TS_ASSERT_EQUALS(daughterCell->GetState<Differentiation>(), Differentiation::TA);

        delete parentCell;
        delete daughterCell; 
    }


    void TestInterModuleCommumication(){

        MAKE_PTR(DummyMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        FixedDurationTestCellCycle* cycleLogic = new FixedDurationTestCellCycle(newCell, CellCycle::M, 0.0, 2.0);
        SimpleAssymmetricStemDivision* diffLogic = new SimpleAssymmetricStemDivision(newCell, Differentiation::S, 0);
        newCell->SetLogic<CellCycle>(cycleLogic);
        newCell->SetLogic<Differentiation>(diffLogic);

        for(int i=0; i<4; i++){
            
            int oldState = newCell->GetState<Differentiation>();

            newCell->UpdateLogic();
            
            int newState = newCell->GetState<Differentiation>(); 

            if(newState == Differentiation::TA && oldState == Differentiation::S){
                //i.e. differentiation occurs on entry into G1
                TS_ASSERT_EQUALS(newCell->GetState<CellCycle>(), CellCycle::G1);
            }
        }

        delete newCell; 
    }


    void TestCellInboundMessages(){

        MAKE_PTR(DummyMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        DifferentiateIfSignalAboveThresh* diffLogic = new DifferentiateIfSignalAboveThresh(newCell, Differentiation::S);
        diffLogic->SetThresh(5);
        newCell->SetLogic<Differentiation>(diffLogic);

        float signalStrength = 0;
        for(int i=0; i<51; i++) {
            signalStrength += 0.1;
            newCell->SendEnvironmentalSignal<DifferentiationDistanceSignal>(signalStrength);
            newCell->UpdateLogic();
        }

        TS_ASSERT_EQUALS(newCell->GetState<Differentiation>(), Differentiation::D);

        delete newCell; 
    }


    void TestCellOutboundMessages(){

        MAKE_PTR(DummyMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        CellCycleTestWithDivisionMechanism* cycleLogic = new CellCycleTestWithDivisionMechanism(newCell, CellCycle::G1, 0.0, 1.0);
        newCell->SetLogic<CellCycle>(cycleLogic);

        for(int i=0; i<10; i++) {
            
            newCell->UpdateLogic();

            std::vector<int> queuedActions = newCell->GetPendingActions();

            if(newCell->GetState<CellCycle>() == CellCycle::M){
                TS_ASSERT_EQUALS(queuedActions[0], CellActions::CallForDivision);
            }else{
                TS_ASSERT_EQUALS(queuedActions.empty(), true);
            }

            newCell->ClearPendingActions();
        }

        delete newCell; 
    }


    void TestTriggerLogicReplacementWarning(){

        MAKE_PTR(DummyMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        FixedDurationTestCellCycle* cycleLogic = new FixedDurationTestCellCycle(newCell, CellCycle::M, 0.0, 2.0);

        newCell->SetLogic<CellCycle>(cycleLogic);
        newCell->SetLogic<CellCycle>(cycleLogic);

        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "There cannot be two logics of type CellCycle. The original logic has been replaced.");
        Warnings::QuietDestroy();
        delete newCell; 
    }
    
};

#endif //TESTLOGICCELL_HPP
