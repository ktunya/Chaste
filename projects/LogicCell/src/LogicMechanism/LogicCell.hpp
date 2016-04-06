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

#ifndef LOGICCELL_HPP
#define LOGICCELL_HPP

#include "Cell.hpp"
#include "AbstractCellLogic.hpp"
#include "PrintableNames.hpp"
#include "LogicTypes.hpp"
#include <vector>
#include <string>
#include <iostream>


class LogicCell : public boost::enable_shared_from_this<LogicCell>, public Cell {

private:

    double mBirthTime;

    std::vector<AbstractCellLogic*> cellLogicModules;
    std::vector<float> environmentalSignals;
    std::vector<int> pendingOutboundActions;

    void setLogicById(int positionInVector, AbstractCellLogic* logic);


    /** Needed for serialization. */
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<Cell>(*this);
        archive & mBirthTime;
        archive & cellLogicModules;
        archive & environmentalSignals;
        archive & pendingOutboundActions;
    }

public:

    LogicCell(boost::shared_ptr<AbstractCellProperty> pDummyMut, bool archiving=false, CellPropertyCollection prop = CellPropertyCollection() );

    virtual ~LogicCell();

    void SetBirthTime(double birthTime);

    double GetBirthTime() const;

    void InitialiseCellCycleModel();

    double GetAge() const;


    void SetCellCycleModel(AbstractCellCycleModel* pCellCycleModel);

    void SetMutationState(boost::shared_ptr<AbstractCellProperty> pMutationState);

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    bool ReadyToDivide();

    CellPtr Divide();

    LogicCell* MakeNewCell();

    void UpdateLogic();

    void dumpStates();

    std::vector<int> GetPendingActions();
    void ProcessPendingActions();
    void ClearPendingActions();
    void SetPendingAction(int actionCode);


    template<typename LOGIC_TYPE>
    void setLogic(AbstractCellLogic* newLogic)
    {
        int logicTypeId = LOGIC_TYPE::id();

        if ((int) cellLogicModules.size() < logicTypeId + 1)
        {
            cellLogicModules.resize(logicTypeId + 1, NULL);
        }

        if (cellLogicModules[logicTypeId] != NULL)
        {
            std::cout << "WARNING: there cannot be two logics of type " << PRINTABLE_NAME(LOGIC_TYPE) <<
            ". The original logic has been replaced." << std::endl;
        }

        newLogic->setCategory(PRINTABLE_NAME(LOGIC_TYPE));
        cellLogicModules[logicTypeId] = newLogic;
    }

    template<typename LOGIC_TYPE>
    int getState()
    {
        int logicTypeId = LOGIC_TYPE::id();
        AbstractCellLogic* selectedLogicModule;

        if (logicTypeId != -1)
        {
            selectedLogicModule = cellLogicModules[logicTypeId];
        }
        else
        {
            throw std::invalid_argument(std::string("This cell has no logic of type ") + PRINTABLE_NAME(LOGIC_TYPE));
        }

        int state = selectedLogicModule->getState();
        return state;
    }


    template<typename SIGNAL_TYPE>
    void sendEnvironmentalSignal(float inboundSignal){

        int signalTypeId = SIGNAL_TYPE::id();

        if ((int) environmentalSignals.size() < signalTypeId + 1)
        {
            environmentalSignals.resize(signalTypeId + 1, 0.0);
        }

        environmentalSignals[signalTypeId] = inboundSignal;
    }


    template<typename SIGNAL_TYPE>
    float getEnvironmentalSignal(){

        int signalTypeId = SIGNAL_TYPE::id();

        if (signalTypeId != -1){
            return environmentalSignals[signalTypeId];
        }

        std::cout << "WARNING: There is no environmental signal of type " << PRINTABLE_NAME(SIGNAL_TYPE) <<
        ". Returning a default signal value of 0." << std::endl;
        return 0;
    }


};

#endif /*LOGICCELL_HPP_*/
