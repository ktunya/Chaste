#ifndef BETTERBACKWARDEULERIVPODESOLVER_HPP_
#define BETTERBACKWARDEULERIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"
#include <cassert>
#include <vector>

template<int SIZE>
class BetterBackwardEulerIvpOdeSolver  : public AbstractOneStepIvpOdeSolver
{
private:
    /** the epsilon to use in calculating numerical jacobians */
    double mNumericalJacobianEpsilon;
    bool mForceUseOfNumericalJacobian;
    
    /** Working memory : residual vector */
    double mResidual[SIZE];
    /** Working memory : Jacobian matrix */
    double** mJacobian;
    /** Working memory : update vector */
    double mUpdate[SIZE];


    void ComputeResidual(AbstractOdeSystem* pAbstractOdeSystem,
                         double timeStep,
                         double time,
                         std::vector<double>& currentYValues,
                         std::vector<double>& currentGuess)
    {
        std::vector<double> dy(SIZE);//For JC to optimize
        pAbstractOdeSystem->EvaluateYDerivatives(time, currentGuess, dy);
        for(unsigned i=0; i<SIZE; i++)
        {
            
            mResidual[i] = currentGuess[i] - timeStep * dy[i] - currentYValues[i];        
        }                     
    }                    
    

    void ComputeJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                         double timeStep,
                         double time,
                         std::vector<double>& currentYValues,
                         std::vector<double>& currentGuess)
    {
        for(unsigned i = 0 ; i<SIZE ; i++)
        {
            for(unsigned j = 0 ; j <SIZE ; j++)
            {
                mJacobian[i][j]=0.0;
            }   
        }

        if(pAbstractOdeSystem->GetUseAnalytic() && !mForceUseOfNumericalJacobian)
        {
            // the ode system has an analytic jacobian, use that
            AbstractOdeSystemWithAnalyticJacobian *p_ode_system 
                 = static_cast<AbstractOdeSystemWithAnalyticJacobian*>(pAbstractOdeSystem);
            p_ode_system->BetterAnalyticJacobian(currentGuess, mJacobian, time, timeStep);
        }
        else
        {
            ComputeNumericalJacobian(pAbstractOdeSystem,
                                     timeStep,
                                     time,
                                     currentYValues,
                                     currentGuess);
        }
    }
    
    void SolveLinearSystem()
    {
        double fact;
        for (unsigned i=0; i<SIZE; i++)
        {
            for (unsigned ii=i+1; ii<SIZE; ii++)
            {
                fact = mJacobian[ii][i]/mJacobian[i][i];
                for (unsigned j=i; j<SIZE; j++)
                {
                    mJacobian[ii][j] -= fact*mJacobian[i][j];
                }
                mResidual[ii] -= fact*mResidual[i];
            }
        }
        for (int i=SIZE-1; i>-1; i--)
        {
            mUpdate[i] = mResidual[i];
            for (unsigned j=i+1; j<SIZE; j++)
            {
                mUpdate[i] -= mJacobian[i][j]*mUpdate[j];
            }
            mUpdate[i] /= mJacobian[i][i];
        }
    }
    
    double ComputeNorm(double vector[SIZE])
    {
        double norm = 0.0;
        for (unsigned i=0; i<SIZE; i++)
        {
            if (fabs(vector[i]) > norm)
            {
                norm = fabs(vector[i]);
            }
        }
        return norm;
    }
    
    
    void ComputeNumericalJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                                  double timeStep,
                                  double time,
                                  std::vector<double>& currentYValues,
                                  std::vector<double>& currentGuess)
    {
        std::vector<double> residual(SIZE);
        std::vector<double> residual_perturbed(SIZE);
        std::vector<double> guess_perturbed(SIZE);
        
        double epsilon= mNumericalJacobianEpsilon;

        ComputeResidual(pAbstractOdeSystem, timeStep, time, currentYValues, currentGuess);
        for(unsigned i=0; i<SIZE; i++)
        {
            residual[i] = mResidual[i];
        }


        for (unsigned global_column=0; global_column<SIZE; global_column++)
        {
            for(unsigned i=0; i<SIZE; i++)
            {
                guess_perturbed[i] = currentGuess[i];
            }
            
            guess_perturbed[global_column] += epsilon;
 
            ComputeResidual(pAbstractOdeSystem, timeStep, time, currentYValues, guess_perturbed);
            for(unsigned i=0; i<SIZE; i++)
            {
                residual_perturbed[i] = mResidual[i];
            }
            
            // compute residual_perturbed - residual
            double one_over_eps=1.0/epsilon;
            for(unsigned i=0; i<SIZE; i++)
            {
                mJacobian[i][global_column] = one_over_eps*(residual_perturbed[i] - residual[i]);
            }
        }
    }

    


protected:
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& currentYValues,
                             std::vector<double>& nextYValues)
    {
        // check the size of the ode system matches the solvers expected
        assert(SIZE == pAbstractOdeSystem->GetNumberOfStateVariables());
        
        
        unsigned counter = 0;
//        const double eps = 1e-6 * rCurrentGuess[0]; // Our tolerance (should use min(guess) perhaps?)
        const double eps = 1e-6; // JonW tolerance
        double norm = 2*eps;
        
        std::vector<double> current_guess(SIZE);
        current_guess.assign(currentYValues.begin(), currentYValues.end());
        
        while (norm > eps)
        {
            // Calculate Jacobian and mResidual for current guess
            ComputeResidual(pAbstractOdeSystem, timeStep, time, currentYValues, current_guess);
            ComputeJacobian(pAbstractOdeSystem, timeStep, time, currentYValues, current_guess);
//            // Update norm (our style)
//            norm = ComputeNorm(mResidual);
            
            // Solve Newton linear system
            SolveLinearSystem();
    
            // Update norm (JonW style)
            norm = ComputeNorm(mUpdate);
            
            // Update current guess
            for (unsigned i=0; i<SIZE; i++)
            {
                current_guess[i] -= mUpdate[i];
            }
                    
            assert(counter++ < 15); // avoid infinite loops
        }
        nextYValues.assign(current_guess.begin(), current_guess.end());
        
    }    

public:


    
    BetterBackwardEulerIvpOdeSolver()
    {
        // default epsilon
        mNumericalJacobianEpsilon = 1e-6;
        mForceUseOfNumericalJacobian = false;
        mJacobian = new double*[SIZE];
        for(unsigned i=0 ; i<SIZE ; i++)
        {
            mJacobian[i] = new double[SIZE];
        }
    }
    
    ~BetterBackwardEulerIvpOdeSolver()
    {
        for(unsigned i=0 ; i<SIZE ; i++)
        {
            delete mJacobian[i];
        }
        delete mJacobian;
    }
    
    void SetEpsilonForNumericalJacobian(double epsilon)
    {
        assert(epsilon > 0);
        mNumericalJacobianEpsilon = epsilon;
    }
     
    /** Force the solver to use the numerical Jacobian even if the 
     *  ode system is one with an analytical jacobian provided
     */
    void ForceUseOfNumericalJacobian()
    {
        mForceUseOfNumericalJacobian = true;
    }    
};


#endif /*BETTERBACKWARDEULERIVPODESOLVER_HPP_*/
