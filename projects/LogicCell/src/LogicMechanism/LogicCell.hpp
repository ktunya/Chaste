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
#include "Warnings.hpp"
#include <vector>
#include <string>


class LogicCell : public boost::enable_shared_from_this<LogicCell>, public Cell {

private:

    // Since LogicCells do not contain a Chaste Cell Cycle Model, the cell birth time
    // must be stored in the cell, and setters and getters overriden to point to this variable
    double mBirthTime;

    // A vector containing the cell's collection of logic modules. Each module deals with a 
    // specific intracellular function, such as the cell cycle or differentiation 
    std::vector<AbstractCellLogic*> cellLogicModules;

    // A vector containing inbound signals from the environment.
    std::vector<float> environmentalSignals;

    // A vector containing pending actions, usually represented by enums. Currently the only
    // cell action implemented is CallForDivision
    std::vector<int> pendingOutboundActions;

    // Allows a daughter cell's vector of logic modules to be set safely. 
    void SetLogicById(int positionInVector, AbstractCellLogic* logic);


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

    // Essentially an "empty" constructor; only requires a pointer to a dummy mutation state to be provided,
    // which is never referenced. Allows a daughter cell to be constructed without any internal logic, with
    // logic then added later based on parent properties.
    LogicCell(boost::shared_ptr<AbstractCellProperty> pDummyMut, bool archiving=false, CellPropertyCollection prop = CellPropertyCollection() );
    virtual ~LogicCell();


    // Overriden birth time setter and getter - refers to the member of this class rather than the cell
    // cycle model. GetAge is overriden for a similar reason.
    virtual void SetBirthTime(double birthTime);
    virtual double GetBirthTime() const;
    virtual double GetAge() const;

    // More functions that are overriden because they're irrelevant for this type of cell. 
    // Some of these now produce a warning to the effect that LogicCells do not use a mutation state
    // or a Chaste Cell Cycle Model 
    virtual void InitialiseCellCycleModel();
    virtual void SetCellCycleModel(AbstractCellCycleModel* pCellCycleModel);
    virtual void SetMutationState(boost::shared_ptr<AbstractCellProperty> pMutationState);

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    // Overriden ReadyToDivide. Doesn't call UpdateCellCyclePhase, but does call UpdateLogic.
    virtual bool ReadyToDivide();
    
    // Overriden Divide. It now calls divide on each of the parent cell's logic modules and passes
    // the results to the daughter. 
    virtual CellPtr Divide();

    // Called by Divide. As well as doing the usual cell division steps, it calls divide on each of
    // the parent's logic modules, and gives the dughter logic to the daughter cell.
    LogicCell* MakeNewCell();

    // Calls update on each logic module 
    void UpdateLogic();

    // Returns the vector of pending cell actions
    std::vector<int> GetPendingActions();

    // For now, this just detects whether the actions vector contains CellActions::CallForDivision
    // If so, mCanDivide is set to true, triggering a cell division 
    void ProcessPendingActions();
    
    // Clears the pending actions vector
    void ClearPendingActions();

    // Adds a pending action; the action code will be an enum for readability
    void SetPendingAction(int actionCode);

    // Sets the logic of type LOGIC_TYPE.
    // Presently, LOGIC_TYPE can take the values CellCycle, Differentiation and Growth
    // If a logic of the same type is set twice in a simulation, a warning is produced.
    // Example usage: pCell->SetLogic<CellCycle>(pStochasticCycle)
    //                pCell->SetLogic<Growth>(pLinearGrowthUpToMaxRadius)
    template<typename LOGIC_TYPE>
    void SetLogic(AbstractCellLogic* newLogic)
    {
        int logicTypeId = LOGIC_TYPE::id();

        if ((int) cellLogicModules.size() < logicTypeId + 1)
        {
            cellLogicModules.resize(logicTypeId + 1, NULL);
        }

        if (cellLogicModules[logicTypeId] != NULL)
        {
            WARNING("There cannot be two logics of type " << PRINTABLE_NAME(LOGIC_TYPE) <<". The original logic has been replaced.");
        }

        newLogic->SetCategory(PRINTABLE_NAME(LOGIC_TYPE));
        cellLogicModules[logicTypeId] = newLogic;
    }

    // Gets the state of a particular logic module. For example:
    // pCell->GetState<CellCycle>()
    // returns the current cell cycle phase 
    template<typename LOGIC_TYPE>
    int GetState()
    {
        int logicTypeId = LOGIC_TYPE::id();
        AbstractCellLogic* selectedLogicModule;

        if (logicTypeId != -1)
        {
            selectedLogicModule = cellLogicModules[logicTypeId];
        }
        else
        {
            EXCEPTION("This cell has no logic of type "<<PRINTABLE_NAME(LOGIC_TYPE));
        }

        int state = selectedLogicModule->GetState();
        return state;
    }

    // Sends an environmental signal to this cell. For example: 
    // pCell->SendEnvironmentalSignal<DistalTipCellSignal>(5.0)
    // would send a signal from the DistalTipCell with strength 5.0.
    // Signal types are declared in EnvironmentalSignalTypes.(h/c)pp 
    template<typename SIGNAL_TYPE>
    void SendEnvironmentalSignal(float inboundSignal){

        int signalTypeId = SIGNAL_TYPE::id();

        if ((int) environmentalSignals.size() < signalTypeId + 1)
        {
            environmentalSignals.resize(signalTypeId + 1, 0.0);
        }

        environmentalSignals[signalTypeId] = inboundSignal;
    }

    // Gets the current value of a particular environmental signal. For example:
    // pCell->GetEnvironmentalSignal<DistalTipCellSignal>()
    // would return the current strength of the Distal Tip Cell signal.
    template<typename SIGNAL_TYPE>
    float GetEnvironmentalSignal(){

        int signalTypeId = SIGNAL_TYPE::id();

        if (signalTypeId != -1){
            return environmentalSignals[signalTypeId];
        }

        WARNING("There is no environmental signal of type "<< PRINTABLE_NAME(SIGNAL_TYPE)<<". Returning a default signal value of 0.");
        return 0;
    }
};

#endif /*LOGICCELL_HPP_*/