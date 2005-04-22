/**
 *Abstract IvpOdeSolver class. Sets up variables and functions for a numerical solution 
 * technique for an intial value ODE problem. 
*/
#ifndef _ABSTRACTIVPODESOLVER_HPP_
#define _ABSTRACTIVPODESOLVER_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <vector>

class AbstractIvpOdeSolver
{
	public:
	virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				              double startTime,
				              double endTime,
				              double timeStep,
				              std::vector<double> initialConditions) = 0;
};

#endif //_ABSTRACTIVPODESOLVER_HPP_
