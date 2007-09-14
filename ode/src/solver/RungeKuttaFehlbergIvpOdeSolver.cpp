/**
 * Concrete RungeKuttaFehlbergIvpOdeSolver class.
 */
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <cmath>
#include <cfloat>
//#include <iostream>
#include <vector>

const double smidge=1e-10;

/*
 * PROTECTED FUNCTIONS
 */
 
 /**
  * This algorithm implements the `rkf45' routine using adaptive time stepping to 
  * solve as fast as possible to a given accuracy.
  * 
  * Note the change in the use of the solutions vector to give results when they are evaluated.
  * 
  * 
  */
void RungeKuttaFehlbergIvpOdeSolver::InternalSolve(OdeSolution& rSolution,
                                                AbstractOdeSystem* pOdeSystem,
                                                std::vector<double>& rYValues,
                                                std::vector<double>& rWorkingMemory,
                                                double startTime,
                                                double endTime,
                                                double maxTimeStep,
                                                double minTimeStep,
                                                double tolerance,
                                                bool outputSolution)
{
    const unsigned number_of_variables = pOdeSystem->GetNumberOfStateVariables();
    mError.reserve(number_of_variables);
    mk1.reserve(number_of_variables);
    mk2.reserve(number_of_variables);
    mk3.reserve(number_of_variables);
    mk4.reserve(number_of_variables);
    mk5.reserve(number_of_variables);
    mk6.reserve(number_of_variables);
    myk2.reserve(number_of_variables);
    myk3.reserve(number_of_variables);
    myk4.reserve(number_of_variables);
    myk5.reserve(number_of_variables);
    myk6.reserve(number_of_variables);
    
    double current_time = startTime;
    double time_step = maxTimeStep;
    bool got_to_end = false;
    bool accurate_enough = false;
    unsigned number_of_time_steps = 0;
    
    if (outputSolution)
    {   // Write out ICs
        rSolution.rGetTimes().push_back(current_time);
        rSolution.rGetSolutions().push_back(rYValues);
    }
    
    // should never get here if this bool has been set to true;
    assert(!mStoppingEventOccured);
    while ( !got_to_end )
    {
        //std::cout << "New timestep\n" << std::flush;
        while (!accurate_enough)
        {
            accurate_enough = true; // assume it is OK until we check and find otherwise

            // Function that calls the appropriate one-step solver
            CalculateNextYValue(pOdeSystem,
                                time_step,
                                current_time,
                                rYValues,
                                rWorkingMemory);
              
            // Find the maximum error in this vector.
            double max_error = -DBL_MAX;
            for (unsigned i=0; i<number_of_variables ; i++)
            {
                if (mError[i] > max_error)
                {
                    max_error = mError[i];   
                }
            }
            
            if (max_error > tolerance)
            {// Reject the step-size and do it again.
                accurate_enough = false;
                //std::cout << "Approximation rejected\n" << std::flush;
            }
            else
            {
                // step forward the time since step has now been made
                current_time = current_time + time_step;    
                //std::cout << "Approximation accepted with time step = "<< time_step << "\n" << std::flush;
                //std::cout << "max_error = " << max_error << " tolerance = " << tolerance << "\n" << std::flush;
                if (outputSolution)
                {   // Write out ICs
                    //std::cout << "In solver Time = " << current_time << " y = " << rWorkingMemory[0] << "\n" << std::flush;
                    rSolution.rGetTimes().push_back(current_time);
                    rSolution.rGetSolutions().push_back(rWorkingMemory);
                    number_of_time_steps++;
                }
            }
            
            // Set a new step size based on the accuracy here
            AdjustStepSize(time_step, max_error, tolerance, maxTimeStep, minTimeStep);
        }
        
        // For the next timestep check the step doesn't go past the end...
        if (current_time + time_step > endTime)
        {   // Allow a smaller timestep for the final step.
            time_step = endTime - current_time;
        }
        
        if ( pOdeSystem->CalculateStoppingEvent(current_time,
                                                rWorkingMemory) == true )
        {
            mStoppingTime = current_time;
            mStoppingEventOccured = true;
        }
        
        if (mStoppingEventOccured || current_time>=endTime)
        {
            got_to_end = true;   
        }
        
        // Approximation accepted - put it in rYValues
        rYValues.assign(rWorkingMemory.begin(), rWorkingMemory.end());
        accurate_enough = false; // for next loop.
        //std::cout << "Finished Time Step\n" << std::flush;
    }
    rSolution.SetNumberOfTimeSteps(number_of_time_steps);
}


/**
 * Solves a system of ODEs using the Runge Kutta Fehlberg Adaptive timestep Initial Value Problem Ordinary Differential Equation Solver
 *
 * To be used in the form:
 *
 * RungeKuttaFehlbergIvpOdeSolver mySolver
 *
 * OdeSolution solution=mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, MaximumTimeStep, any_old_double);
 *
 * See documentation for AbstractIvpOdeSolver::Solve()
 * 
 * @return error the difference between the 4th and 5th order approximations.
 */
void RungeKuttaFehlbergIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                  double timeStep,
                                                  double time,
                                                  std::vector<double>& currentYValues,
                                                  std::vector<double>& nextYValues)
{
    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();


    std::vector<double>& dy = nextYValues; // re-use memory (not that it makes much difference here!)
    
    pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValues, dy);
    
    for (unsigned i=0;i<num_equations; i++)
    {
        mk1[i] = timeStep*dy[i];
        myk2[i] = currentYValues[i] + 0.25*mk1[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time + 0.25*timeStep, myk2, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        mk2[i] = timeStep*dy[i];
        myk3[i] = currentYValues[i] + 0.09375*mk1[i] + 0.28125*mk2[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time + 0.375*timeStep, myk3, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        mk3[i] = timeStep*dy[i];
        myk4[i] = currentYValues[i] + m1932o2197*mk1[i] - m7200o2197*mk2[i] 
                    + m7296o2197*mk3[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+m12o13*timeStep, myk4, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        mk4[i] = timeStep*dy[i];
        myk5[i] = currentYValues[i] + m439o216*mk1[i] - 8*mk2[i] 
                    + m3680o513*mk3[i]- m845o4104*mk4[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep, myk5, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        mk5[i] = timeStep*dy[i];
        myk6[i] = currentYValues[i] - m8o27*mk1[i] + 2*mk2[i] - m3544o2565*mk3[i]
                        + m1859o4104*mk4[i] - 0.275*mk5[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep, myk6, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        mk6[i]=timeStep*dy[i];
        mError[i] = (1/timeStep)*fabs(m1o360*mk1[i] - m128o4275*mk3[i]
                    - m2197o75240*mk4[i] + 0.02*mk5[i]+ m2o55*mk6[i]);
        nextYValues[i] = currentYValues[i] + m25o216*mk1[i] + m1408o2565*mk3[i]
                        + m2197o4104*mk4[i] - 0.2*mk5[i];
    }
}

void RungeKuttaFehlbergIvpOdeSolver::AdjustStepSize(double& rCurrentStepSize,
                                const double& rError, 
                                const double& rTolerance, 
                                const double& rMaxTimeStep, 
                                const double& rMinTimeStep)
{
    // Work out scaling factor delta for the step size
    double delta = pow(rTolerance/(2.0*rError), 0.25);
    
    // Maximum adjustment is *0.1 or *4
    if (delta <= 0.1)
    {
        rCurrentStepSize *= 0.1;
    }
    else if (delta >= 4.0)
    {
        rCurrentStepSize *= 4.0;  
    }
    else
    {
        rCurrentStepSize *= delta;
    }
    
    if (rCurrentStepSize > rMaxTimeStep)
    {
        rCurrentStepSize = rMaxTimeStep;
    }
    
    if (rCurrentStepSize < rMinTimeStep)
    {
        EXCEPTION("RKF45 Solver: Ode needs a smaller timestep than the set minimum\n");   
    }
    
}     


/*
 * PUBLIC FUNCTIONS
 */

OdeSolution RungeKuttaFehlbergIvpOdeSolver::Solve(AbstractOdeSystem* pOdeSystem,
                                               std::vector<double>& rYValues,
                                               double startTime,
                                               double endTime,
                                               double timeStep,
                                               double tolerance)
{
    assert(rYValues.size()==pOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);
    
    mStoppingEventOccured = false;
    if ( pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("(Solve with sampling) Stopping event is true for initial condition");
    }
    // Allocate working memory
    std::vector<double> working_memory(rYValues.size());
    // And solve...
    OdeSolution solutions;
    //solutions.SetNumberOfTimeSteps((unsigned)(10.0*(startTime-endTime)/timeStep));
    bool return_solution = true;
    InternalSolve(solutions, pOdeSystem, rYValues, working_memory, startTime, endTime, timeStep, 1e-4, tolerance, return_solution);
    return solutions;
}

void RungeKuttaFehlbergIvpOdeSolver::Solve(AbstractOdeSystem* pOdeSystem,
                                        std::vector<double>& rYValues,
                                        double startTime,
                                        double endTime,
                                        double timeStep)
{
    assert(rYValues.size()==pOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);
    
    mStoppingEventOccured = false;
    if ( pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("(Solve without sampling) Stopping event is true for initial condition");
    }
    // Allocate working memory
    std::vector<double> working_memory(rYValues.size());
    // And solve...
    OdeSolution not_required_solution;
    bool return_solution = false;
    InternalSolve(not_required_solution, pOdeSystem, rYValues, working_memory, startTime, endTime, timeStep, 1e-4, 1e-5, return_solution);
}


