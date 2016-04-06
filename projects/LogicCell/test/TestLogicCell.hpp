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

#include "LogicCell.hpp"
#include "AbstractCellLogic.hpp"
#include "LogicTypes.hpp"
#include "EnvironmentalSignalTypes.hpp"
#include "CellActions.hpp"

#include "FixedDurationCellCycle.hpp"
#include "SimpleAssymmetricStemDivision.hpp"
#include "DifferentiateIfSignalBelowThresh.hpp"
#include "CellCycleWithDivisionMechanism.hpp"

//TODO: Temporary use of this mutation state and cell cycle model. Write a dummy cell cycle model and
//      mutation state instead that throw errors when data is accessed.
#include "WildTypeCellMutationState.hpp"


#ifndef TESTLOGICCELL_HPP
#define TESTLOGICCELL_HPP

class TestLogicCell : public AbstractCellBasedTestSuite {


public:

    void TestInitialisation(){

        MAKE_PTR(WildTypeCellMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        FixedDurationCellCycle* cycleLogic = new FixedDurationCellCycle(newCell, CellCycle::G1, 0.0, 2.0);
        newCell->setLogic<CellCycle>(cycleLogic);

        int state = newCell->getState<CellCycle>();
        TS_ASSERT_EQUALS(state, CellCycle::G1);

        for(int i=0; i<5; i++){
            newCell->update();
        }

        state = newCell->getState<CellCycle>();
        TS_ASSERT_EQUALS(state, CellCycle::G2);

        delete newCell;
    }


    void TestSymmetricDivision(){

        MAKE_PTR(WildTypeCellMutationState, pDummyMut);
        LogicCell* parentCell = new LogicCell(pDummyMut);

        FixedDurationCellCycle* cycleLogic = new FixedDurationCellCycle(parentCell, CellCycle::G1, 0.0, 2.0);
        parentCell->setLogic<CellCycle>(cycleLogic);

        LogicCell* daughterCell = parentCell->MakeNewCell();

        TS_ASSERT_EQUALS(parentCell->getState<CellCycle>(), daughterCell->getState<CellCycle>());

        delete parentCell;
        delete daughterCell;
    }


    void TestAssymmetricDivision(){

        MAKE_PTR(WildTypeCellMutationState, pDummyMut);
        LogicCell* parentCell = new LogicCell(pDummyMut);

        FixedDurationCellCycle* cycleLogic = new FixedDurationCellCycle(parentCell, CellCycle::M, 0.0, 2.0);
        parentCell->setLogic<CellCycle>(cycleLogic);
        SimpleAssymmetricStemDivision* diffLogic = new SimpleAssymmetricStemDivision(parentCell, Differentiation::S);
        parentCell->setLogic<Differentiation>(diffLogic);

        LogicCell* daughterCell = parentCell->MakeNewCell();

        TS_ASSERT_EQUALS(daughterCell->getState<Differentiation>(), Differentiation::TA);

        delete parentCell;
        delete daughterCell;
    }


    void TestInterModuleCommumication(){

        MAKE_PTR(WildTypeCellMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        FixedDurationCellCycle* cycleLogic = new FixedDurationCellCycle(newCell, CellCycle::M, 0.0, 2.0);
        SimpleAssymmetricStemDivision* diffLogic = new SimpleAssymmetricStemDivision(newCell, Differentiation::S);
        newCell->setLogic<CellCycle>(cycleLogic);
        newCell->setLogic<Differentiation>(diffLogic);

        for(int i=0; i<4; i++){
            bool wasStem = false;
            if(newCell->getState<Differentiation>() == Differentiation::S){
                wasStem = true;
            }
            newCell->update();
            if(newCell->getState<Differentiation>() == Differentiation::TA && wasStem){
                TS_ASSERT_EQUALS(newCell->getState<CellCycle>(), CellCycle::G1);
            }
        }

        delete newCell;
    }


    void TestCellInboundMessages(){

        MAKE_PTR(WildTypeCellMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        DifferentiateIfSignalBelowThresh* diffLogic = new DifferentiateIfSignalBelowThresh(newCell, Differentiation::S);
        diffLogic->setThresh(5);

        newCell->setLogic<Differentiation>(diffLogic);

        float signalStrength = 10;
        for(int i=0; i<52; i++) {
            signalStrength -= 0.1;
            newCell->sendEnvironmentalSignal<MyDiffusableSignal>(signalStrength);
            newCell->update();
        }

        TS_ASSERT_EQUALS(newCell->getState<Differentiation>(), Differentiation::D);

        delete newCell;
    }


    void TestCellOutboundMessages(){

        MAKE_PTR(WildTypeCellMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        CellCycleWithDivisionMechanism* cycleLogic = new CellCycleWithDivisionMechanism(newCell, CellCycle::G1, 0.0, 1.0);
        newCell->setLogic<CellCycle>(cycleLogic);

        for(int i=0; i<10; i++) {
            newCell->update();

            std::vector<int> queuedActions = newCell->getPendingActions();

            if(newCell->getState<CellCycle>()==CellCycle::M){
                TS_ASSERT_EQUALS(queuedActions[0], CellActions::CallForDivision);
            }else{
                TS_ASSERT_EQUALS(queuedActions.empty(), true);
            }

            newCell->clearPendingActions();
        }

        delete newCell;
    }


    void TestTriggerLogicReplacementWarning(){

        MAKE_PTR(WildTypeCellMutationState, pDummyMut);
        LogicCell* newCell = new LogicCell(pDummyMut);

        FixedDurationCellCycle* cycleLogic = new FixedDurationCellCycle(newCell, CellCycle::M, 0.0, 2.0);

        newCell->setLogic<CellCycle>(cycleLogic);
        newCell->setLogic<CellCycle>(cycleLogic);

        delete newCell;
    }

};

#endif //TESTLOGICCELL_HPP
