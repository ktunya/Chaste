/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "CellCycleModelOdeHandler.hpp"

CellCycleModelOdeHandler::CellCycleModelOdeHandler(double lastTime,
                                                   boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : mDt(DOUBLE_UNSET),
      mpOdeSystem(NULL),
      mpOdeSolver(pOdeSolver),
      mLastTime(lastTime)
{
}

CellCycleModelOdeHandler::~CellCycleModelOdeHandler()
{
    if (mpOdeSystem != NULL)
    {
        delete mpOdeSystem;
    }
}

void CellCycleModelOdeHandler::SetOdeSystem(AbstractOdeSystem* pOdeSystem)
{
    mpOdeSystem = pOdeSystem;
}

AbstractOdeSystem* CellCycleModelOdeHandler::GetOdeSystem() const
{
    return mpOdeSystem;
}

const boost::shared_ptr<AbstractCellCycleModelOdeSolver> CellCycleModelOdeHandler::GetOdeSolver() const
{
    return mpOdeSolver;
}

void CellCycleModelOdeHandler::SetTimeStep(double timeStep)
{
    mDt = timeStep;
}

double CellCycleModelOdeHandler::GetDt()
{
    if (mDt == DOUBLE_UNSET)
    {
#ifdef CHASTE_CVODE
        mDt = SimulationTime::Instance()->GetTimeStep();
#else
        mDt = 0.0001; // Some models need this, so let's pick a safe default
#endif // CHASTE_CVODE
    }
    return mDt;
}

bool CellCycleModelOdeHandler::SolveOdeToTime(double currentTime)
{
    AdjustOdeParameters(currentTime);

    mpOdeSolver->SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, currentTime, GetDt());

    bool stopping_event_occurred = mpOdeSolver->StoppingEventOccurred();
    if (stopping_event_occurred)
    {
        mLastTime = mpOdeSolver->GetStoppingTime();
    }
    else
    {
        mLastTime = currentTime;
    }
    
    return stopping_event_occurred;
}

void CellCycleModelOdeHandler::AdjustOdeParameters(double currentTime)
{
}
