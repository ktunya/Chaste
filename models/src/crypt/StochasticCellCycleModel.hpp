#ifndef STOCHASTICCELLCYCLEMODEL_HPP_
#define STOCHASTICCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 *  Stochastic cell model
 *
 *  Cell cycle time is deterministic for stem cells and stochastic (normally
 *  distributed with mean CancerParameters::TransitCellCycleTime and variance 1)
 */
class StochasticCellCycleModel : public AbstractCellCycleModel
{
private:
    RandomNumberGenerator* mpGen;
    
public:
    StochasticCellCycleModel(RandomNumberGenerator *pGen)
    {
        mpGen=pGen;
    }
    
    bool ReadyToDivide(double timeSinceBirth);
    
    AbstractCellCycleModel *CreateCellCycleModel();
};

#endif /*STOCHASTICCELLCYCLEMODEL_HPP_*/
