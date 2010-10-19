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

#ifdef CHASTE_CVODE

#include <sstream>
#include <iostream>
#include <cmath>

#include "AbstractCvodeCell.hpp"
#include "CvodeAdaptor.hpp" // For CvodeErrorHandler
#include "TimeStepper.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"

// CVODE headers
#include <cvode/cvode.h>
#include <sundials/sundials_nvector.h>
#include <cvode/cvode_dense.h>

/**
 * Callback function provided to CVODE to allow it to 'call' C++ member functions
 * (in particular, AbstractCvodeCell::EvaluateRhs).
 *
 * @param t  current time
 * @param y  state variable vector
 * @param ydot  derivatives vector to be filled in
 * @param pData  pointer to the cell being simulated
 */
int AbstractCvodeCellRhsAdaptor(realtype t, N_Vector y, N_Vector ydot, void *pData)
{
    assert(pData != NULL);
    AbstractCvodeCell* pCell = (AbstractCvodeCell*) pData;
    try
    {
        pCell->EvaluateRhs(t, y, ydot);
    }
    catch (Exception &e)
    {
        std::cerr << "CVODE RHS Exception: " << e.GetMessage()
                  << std::endl << std::flush;
        return -1;
    }
    return 0;
}


AbstractCvodeCell::AbstractCvodeCell(boost::shared_ptr<AbstractIvpOdeSolver> pSolver /* unused */,
                                     unsigned numberOfStateVariables,
                                     unsigned voltageIndex,
                                     boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCellInterface(pSolver, voltageIndex, pIntracellularStimulus),
      AbstractParameterisedSystem<N_Vector>(numberOfStateVariables),
      mpCvodeMem(NULL),
      mMaxSteps(0)
{
    SetTolerances();
}


AbstractCvodeCell::~AbstractCvodeCell()
{
    if (mStateVariables)
    {
        mStateVariables->ops->nvdestroy(mStateVariables);
        mStateVariables = NULL;
    }
    if (mParameters)
    {
        mParameters->ops->nvdestroy(mParameters);
        mParameters = NULL;
    }
}


void AbstractCvodeCell::Init()
{
    SetStateVariables(GetInitialConditions());
    if (mParameters)
    {
        mParameters->ops->nvdestroy(mParameters);
        mParameters = NULL;
    }
    mParameters = N_VNew_Serial(rGetParameterNames().size());
    for (unsigned i=0; i<NV_LENGTH_S(mParameters); i++)
    {
        NV_Ith_S(mParameters, i) = 0.0;
    }
}


double AbstractCvodeCell::GetVoltage()
{
    assert(mStateVariables);
    return NV_Ith_S(mStateVariables, mVoltageIndex);
}

void AbstractCvodeCell::SetVoltage(double voltage)
{
    assert(mStateVariables);
    NV_Ith_S(mStateVariables, mVoltageIndex) = voltage;
}


void AbstractCvodeCell::ResetToInitialConditions()
{
    SetStateVariables(GetInitialConditions());
}

N_Vector AbstractCvodeCell::GetInitialConditions()
{
    assert(mpSystemInfo);
    std::vector<double> inits = mpSystemInfo->GetInitialConditions();
    N_Vector v = N_VNew_Serial(inits.size());
    for (unsigned i=0; i<inits.size(); i++)
    {
        NV_Ith_S(v, i) = inits[i];
    }
    return v;
}

N_Vector AbstractCvodeCell::GetStateVariables()
{
    assert(mpSystemInfo);
    assert(mStateVariables);
    return CopyVector(mStateVariables);
}

void AbstractCvodeCell::SetStateVariables(N_Vector stateVars)
{
    if (mStateVariables and stateVars != mStateVariables)
    {
        mStateVariables->ops->nvdestroy(mStateVariables);
        ///\todo #890 re-init CVODE here?
    }
    mStateVariables = stateVars;
}

void AbstractCvodeCell::SetStateVariablesUsingACopyOfThisVector(N_Vector stateVars)
{
    // Cope if stateVars == mStateVariables
    N_Vector temp = CopyVector(stateVars);
    if (mStateVariables)
    {
        mStateVariables->ops->nvdestroy(mStateVariables);
    }
    mStateVariables = temp;
}


void AbstractCvodeCell::SetVoltageDerivativeToZero(bool clamp)
{
    mSetVoltageDerivativeToZero = clamp;
}


OdeSolution AbstractCvodeCell::Compute(double tStart, double tEnd, double tSamp)
{
    if (tSamp == 0.0)
    {
        tSamp = HeartConfig::Instance()->GetPrintingTimeStep();
    }
    return Solve(tStart, tEnd, tSamp, tSamp);
}


void AbstractCvodeCell::ComputeExceptVoltage(double tStart, double tEnd)
{
    EXCEPTION("This method is not yet implemented for CVODE cells.");
}


OdeSolution AbstractCvodeCell::Solve(realtype tStart,
                                     realtype tEnd,
                                     realtype maxDt,
                                     realtype tSamp)
{
    assert(tEnd >= tStart);
    assert(tSamp > 0.0);

    SetupCvode(mStateVariables, tStart, maxDt);

    TimeStepper stepper(tStart, tEnd, tSamp);

    // Set up ODE solution
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(MakeStdVec(mStateVariables));
    solutions.rGetTimes().push_back(tStart);
    solutions.SetOdeSystemInformation(mpSystemInfo);

    // Main time sampling loop
    while (!stepper.IsTimeAtEnd())
    {
//        std::cout << "Solving to time " << stepper.GetNextTime() << std::endl << std::flush;
        double cvode_stopped_at;
        int ierr;
        try
        {
            ierr = CVode(mpCvodeMem, stepper.GetNextTime(), mStateVariables,
                         &cvode_stopped_at, CV_NORMAL);
        }
        catch (...)
        {
            FreeCvodeMemory();
            throw;
        }
        if (ierr<0) CvodeError(ierr, "CVODE failed to solve system");
        // Not root finding, so should have reached requested time
        assert(fabs(cvode_stopped_at - stepper.GetNextTime()) < DBL_EPSILON);

        VerifyStateVariables();

        // Store solution
        solutions.rGetSolutions().push_back(MakeStdVec(mStateVariables));
        solutions.rGetTimes().push_back(cvode_stopped_at);
        stepper.AdvanceOneTimeStep();
    }

    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTotalTimeStepsTaken());

    int ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning
    FreeCvodeMemory();

    return solutions;
}

void AbstractCvodeCell::Solve(realtype tStart,
                              realtype tEnd,
                              realtype maxDt)
{
    assert(tEnd >= tStart);

    SetupCvode(mStateVariables, tStart, maxDt);

    double cvode_stopped_at;
    int ierr;
    try
    {
        ierr = CVode(mpCvodeMem, tEnd, mStateVariables, &cvode_stopped_at, CV_NORMAL);
    }
    catch (...)
    {
        FreeCvodeMemory();
        throw;
    }
    if (ierr<0) CvodeError(ierr, "CVODE failed to solve system");
    // Not root finding, so should have reached requested time
    assert(fabs(cvode_stopped_at - tEnd) < DBL_EPSILON);

    ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning
    FreeCvodeMemory();

    VerifyStateVariables();
}


void AbstractCvodeCell::SetMaxSteps(long int numSteps)
{
    mMaxSteps = numSteps;
}

long int AbstractCvodeCell::GetMaxSteps()
{
    return mMaxSteps;
}

void AbstractCvodeCell::SetTolerances(double relTol, double absTol)
{
    mRelTol = relTol;
    mAbsTol = absTol;
}

double AbstractCvodeCell::GetRelativeTolerance()
{
    return mRelTol;
}

double AbstractCvodeCell::GetAbsoluteTolerance()
{
    return mAbsTol;
}

double AbstractCvodeCell::GetLastStepSize()
{
    return mLastInternalStepSize;
}


void AbstractCvodeCell::SetupCvode(N_Vector initialConditions,
                                   realtype tStart,
                                   realtype maxDt)
{
    assert(NV_LENGTH_S(initialConditions) == GetNumberOfStateVariables());
    assert(maxDt >= 0.0);

    mpCvodeMem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (mpCvodeMem == NULL) EXCEPTION("Failed to SetupCvode CVODE");
    // Set error handler
    CVodeSetErrHandlerFn(mpCvodeMem, CvodeErrorHandler, NULL);
    // Set the user data
    CVodeSetFdata(mpCvodeMem, (void*)(this));
    // Setup CVODE
    CVodeMalloc(mpCvodeMem, AbstractCvodeCellRhsAdaptor, tStart, initialConditions,
                CV_SS, mRelTol, &mAbsTol);
    CVodeSetMaxStep(mpCvodeMem, maxDt);
    // Change max steps if wanted
    if (mMaxSteps > 0)
    {
        CVodeSetMaxNumSteps(mpCvodeMem, mMaxSteps);
    }
    // Attach a linear solver for Newton iteration
    CVDense(mpCvodeMem, NV_LENGTH_S(initialConditions));
}


void AbstractCvodeCell::FreeCvodeMemory()
{
    CVodeFree(&mpCvodeMem);
}


// Errors get caught before they hit this, but we keep it just in case...
#define COVERAGE_IGNORE
void AbstractCvodeCell::CvodeError(int flag, const char * msg)
{
    std::stringstream err;
    err << msg << ": " << CVodeGetReturnFlagName(flag);
    std::cerr << err.str() << std::endl << std::flush;
    EXCEPTION(err.str());
}
#undef COVERAGE_IGNORE


std::vector<double> AbstractCvodeCell::MakeStdVec(N_Vector v)
{
    unsigned size = NV_LENGTH_S(v);
    std::vector<double> sv(size);
    for (unsigned i=0; i<size; i++)
    {
        sv[i] = NV_Ith_S(v, i);
    }
    return sv;
}


std::string AbstractCvodeCell::DumpState(const std::string& message,
                                         N_Vector Y)
{
    std::stringstream res;
    res << message << "\nState:\n";
    if (Y == NULL)
    {
        Y = mStateVariables;
    }
    for (unsigned i=0; i<NV_LENGTH_S(Y); i++)
    {
        res << "\t" << rGetStateVariableNames()[i] << ":" << NV_Ith_S(Y, i) << "\n";
    }
    return res.str();
}

N_Vector AbstractCvodeCell::CopyVector(N_Vector originalVec)
{
    unsigned size = NV_LENGTH_S(originalVec);
    N_Vector v = N_VClone(originalVec);
    for (unsigned i=0; i<size; i++)
    {
        NV_Ith_S(v, i) = NV_Ith_S(originalVec, i);
    }
    return v;
}


#endif // CHASTE_CVODE
