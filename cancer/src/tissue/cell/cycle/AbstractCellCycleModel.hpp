#ifndef ABSTRACTCELLCYCLEMODEL_HPP_
#define ABSTRACTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "CellTypes.hpp"
#include "CellCyclePhases.hpp"
#include "SimulationTime.hpp"
#include "CancerParameters.hpp"
#include "TissueCell.hpp"
#include <vector>

// Needs to be included last
#include <boost/serialization/export.hpp>

class TissueCell; // Circular definition (cells need to know about cycle models and vice-versa).

/**
 * The AbstractCellCycleModel contains basic information to all cell cycle models.
 * It handles assignment of birth time, cell cycle phase and a TissueCell.
 */
class AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mBirthTime;

        // Make sure the simulation and cancer parameters get saved too
        SimulationTime* p_time = SimulationTime::Instance();
        archive & *p_time;
        archive & p_time;
        CancerParameters* p_params = CancerParameters::Instance();
        archive & *p_params;
        archive & p_params;

        archive & mGeneration;
        // DO NOT archive & mpCell; -- The CellCycleModel is only ever archived from the Cell 
        // which knows this and it is handled in the load_construct of TissueCell.
        archive & mCurrentCellCyclePhase;
    }

        
protected:
    /** The cell that this model is associated with */
    TissueCell* mpCell;
    
    /** 
     * The time that the cell began to split from its parent
     * (i.e. beginning of M phase NOT the end) 
     */
    double mBirthTime;
    
    /** The phase of the cell cycle that this model is in (specified in CellCyclePhases.hpp) */
    CellCyclePhase mCurrentCellCyclePhase;
    
    /** The Generation of this cell (STEM cells have a generation of 0) */
    unsigned mGeneration;    
        
public:

    /**
     * Sets up a new AbstractCellCycleModel, gives it a birth time of the 
     * current simulation time (which is overwritten by some subclasses) 
     */
    AbstractCellCycleModel()
        : mpCell(NULL),
          mBirthTime(SimulationTime::Instance()->GetDimensionalisedTime()),
          mCurrentCellCyclePhase(M_PHASE) 
    {
        mGeneration = 0;
    }

    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractCellCycleModel();
    
    virtual void SetCell(TissueCell* pCell);
    
    virtual void Initialise() {};
    
    TissueCell* GetCell();
    
    /**
     * Set the cell's time of birth (usually not required as it should be inside
     * the indivdual cell-cycle-model-constructor, but useful for tests).
     * 
     * @param birthTime the simulation time at this cell's birth.
     * 
     * (This function is overridden in AbstractOdeBasedCellCycleModel).
     */
    virtual void SetBirthTime(double birthTime);
    
    /**
     * @return the time at which the cell was born.
     */
    double GetBirthTime() const;
    
    /**
     * Returns the cell's age...
     */
    double GetAge();
    
    /**
     * Sets the cell's generation...
     */
    void SetGeneration(unsigned generation);    
    
    /**
     * Returns the cell's generation...
     */
    unsigned GetGeneration() const;

    
    /**
     * Returns the cell types of the next generation of cells in a vector
     * [0] is the new mother cell type, [1] is the new daughter cell type
     *
     * \todo perhaps returning a std::pair would be better?
     */
    virtual std::vector<CellType> GetNewCellTypes();
    
    /**
     * Determine whether the cell is ready to divide (enter M phase).
     *
     * The intention is that this method is called precisely once at
     * each timestep of the simulation.  However this does not appear
     * to always be the case at present.
     */
    virtual bool ReadyToDivide()=0;
    
    /**
     * Each cell cycle model must be able to be reset 'after' a cell division.
     * 
     * Actually, this method is called from TissueCell::Divide to
     * reset the cell cycle just before the daughter cell is created.
     * CreateCellCycleModel can then clone our state to generate a
     * cell cycle model instance for the daughter cell.
     */
    virtual void ResetModel()=0;
    
    /**
     * Builder method to create new instances of the cell cycle model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     * 
     * This method is called in 2 circumstances:
     *  - By the copy constructor and operator= of TissueCell to create a copy of the cell cycle
     *    model when copying a cell.  This method therefore just needs to create any instance
     *    of the right class, as operator= on the cell cycle model is then called to ensure the
     *    model is copied properly.
     *  - By TissueCell.Divide to create a cell cycle model for the daughter cell.
     *    This method must thus produce a cell cycle model in a suitable state for a
     *    newly-born cell spawned from the 'current' cell.
     *    Note that the parent cell cycle model will have had ResetModel() called just
     *    before CreateCellCycleModel is called.
     */
    virtual AbstractCellCycleModel *CreateCellCycleModel()=0;
    
    /**
     * This normally does nothing but is over-ridden when the mother cell has
     * an AbstractSimpleMeineke cell cycle model and cell is a stem.
     */ 
    virtual void SetMotherGeneration();
    

    /**
     * @return whether the cell cycle model uses beta-catenin levels in cell cycle model, ie Inge models.
     */
    virtual bool UsesBetaCat();

    /**
     * @return the level of membrane bound beta-catenin. However in most Cell Cycle models this does not exist.
     * We have a "work-around" such that we throw an error if we try and access it for any other cell type.
     * \todo may be better to use dynamic_cast and/or MI.
     */
    virtual double GetMembraneBoundBetaCateninLevel();
    
    /**
     * @return the level of cytoplasm beta-catenin. However in most Cell Cycle models this does not exist.
     * We have a "work-around" such that we throw an error if we try and access it for any other cell type. 
     * \todo may be better to use dynamic_cast and/or MI.
     */
    virtual double GetCytoplasmicBetaCateninLevel();
    
    /**
     * @return the level of nuclear bound beta-catenin. However in most Cell Cycle models this does not exist.
     * We have a "work-around" such that we throw an error if we try and access it for any other cell type. 
     * \todo may be better to use dynamic_cast and/or MI.
     */
    virtual double GetNuclearBetaCateninLevel();
    

    /*
     * @return the current cell cycle phase
     */
    CellCyclePhase GetCurrentCellCyclePhase();
    
    /**
     * @return the duration of the S phase of the cell cycle
     */
    virtual double GetSDuration();
    
    /**
     * @return the duration of the G2 phase of the cell cycle
     */
    virtual double GetG2Duration();
    
    /**
     * @return the duration of the M phase of the cell cycle
     */
    virtual double GetMDuration();
    
};


BOOST_IS_ABSTRACT(AbstractCellCycleModel)


#endif /*ABSTRACTCELLCYCLEMODEL_HPP_*/
