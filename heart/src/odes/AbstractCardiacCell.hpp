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

#ifndef ABSTRACTCARDIACCELL_HPP_
#define ABSTRACTCARDIACCELL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/version.hpp>

#include "AbstractOdeSystem.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"

#include <vector>

typedef enum _CellModelState
{
    STATE_UNSET = 0,
    FAST_VARS_ONLY,
    ALL_VARS
} CellModelState;

/**
 * This is the base class for ode-based cardiac cell models.
 *
 * It is essentially a cardiac-specific wrapper for ODE systems
 * providing an interface which can interact with the Stimulus
 * classes and the voltage in a mono/bidomain simulation.
 *
 * Concrete classes can be autogenerated from CellML files
 * by the pyCML package, and will automatically inherit from this class.
 */
class AbstractCardiacCell : public AbstractOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
        archive & mDt;
        archive & mSetVoltageDerivativeToZero;
        if (version > 0)
        {
            // Note that when loading a version 0 archive, this will be initialised to
            // false by our constructor.  So we should get a consistent (wrong) answer
            // with previous versions of Chaste when in tissue.
            archive & mIsUsedInTissue;
        }
        // archive & mVoltageIndex; - always set by constructor - called by concrete class
        // archive & mpOdeSolver; - always set by constructor - called by concrete class
        // archive & mpIntracellularStimulus; - always set by constructor - called by concrete class
    }

protected:
    /** The index of the voltage within our state variable vector. */
    unsigned mVoltageIndex;
    /** Pointer to the solver used to simulate this cell. */
    boost::shared_ptr<AbstractIvpOdeSolver> mpOdeSolver;
    /** The timestep to use when simulating this cell.  Set from the HeartConfig object. */
    double mDt;
    /** The intracellular stimulus current. */
    boost::shared_ptr<AbstractStimulusFunction> mpIntracellularStimulus;

    /**
     * Flag set to true if ComputeExceptVoltage is called, to indicate
     * to subclass EvaluateYDerivatives methods that V should be
     * considered fixed, and hence dV/dt set to zero.
     */
    bool mSetVoltageDerivativeToZero;

    /** Whether this cell exists in a tissue, or is an isolated cell. */
    bool mIsUsedInTissue;

public:
    /** Create a new cardiac cell. The state variables of the cell will be 
     *  set to AbstractOdeSystemInformation::GetInitialConditions(). Note that
     *  calls to SetDefaultInitialConditions() on a particular instance of this class
     *  will not modify its state variables. You can modify them directly with 
     *  rGetStateVariables(). 
     *
     * @param pOdeSolver  the ODE solver to use when simulating this cell
     * @param numberOfStateVariables  the size of the ODE system modelling this cell
     * @param voltageIndex  the index of the transmembrane potential within the vector of state variables
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    AbstractCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                        unsigned numberOfStateVariables,
                        unsigned voltageIndex,
                        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /** Virtual destructor */
    virtual ~AbstractCardiacCell();

    /**
     * Initialise the cell:
     *  - set our state variables to the initial conditions,
     *  - resize model parameters vector.
     * 
     * @note Must be called by subclass constructors.
     */
    void Init();

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual OdeSolution Compute(double tStart, double tEnd);

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * with timestep #mDt, but does not update the voltage.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void ComputeExceptVoltage(double tStart, double tEnd);

    /**
     * Computes the total current flowing through the cell membrane, using the current
     * values of the state variables.
     *
     * Should return a value in units of microamps/cm^2.  Note that many cell models
     * do not use these dimensions (let alone these units) and so a complex conversion
     * is required.  There are 2 main cases:
     *   - Cell model uses amps per unit capacitance.  Often in this case the units used
     *     for the cell capacitance don't make sense (e.g. uF/cm^2 is used, and dV/dt=I/C_m).
     *     Hence we suggest examining the equation for dV/dt given in the model to determine
     *     what the cell model really considers the value for C_m to be, and scaling by
     *     Chaste's C_m / cell model C_m (the latter implicitly being dimensionless).
     *   - Cell model uses amps.  In this case you need to divide by an estimate of the cell
     *     surface area.  Assuming the model represents a single cell, and gives C_m in farads,
     *     then scaling by Chaste's C_m / model C_m seems reasonable.  If the 'cell model'
     *     doesn't actually represent a single whole cell, then you'll have to be more careful.
     * In both cases additional scaling may be required to obtain correct units once the
     * dimensions have been sorted out.
     *
     * Chaste's value for C_m can be obtained from HeartConfig::Instance()->GetCapacitance()
     * and is measured in uF/cm^2.
     */
    virtual double GetIIonic() = 0;

    /** Set the transmembrane potential
     * @param voltage  new value
     */
    void SetVoltage(double voltage);

    /**
     * Get the current value of the transmembrane potential, as given
     * in our state variable vector.
     */
    double GetVoltage();

    /** Get the index of the transmembrane potential within our state variable vector. */
    unsigned GetVoltageIndex();

    /**
     * Set the intracellular stimulus.
     * Shorthand for SetIntracellularStimulusFunction.
     * @param pStimulus  new stimulus function
     */
    void SetStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus);

    /**
     * Get the value of the intracellular stimulus.
     * Shorthand for GetIntracellularStimulus.
     * @param time  the time at which to evaluate the stimulus
     */
    double GetStimulus(double time);

    /**
     * Set the intracellular stimulus.
     * This should have units of uA/cm^2 for single-cell problems,
     * or uA/cm^3 in a tissue simulation.
     * @param pStimulus  new stimulus function
     */
    void SetIntracellularStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus);

    /**
     * Get the value of the intracellular stimulus.
     * This will have units of uA/cm^2 for single-cell problems,
     * or uA/cm^3 in a tissue simulation.
     *
     * @param time  the time at which to evaluate the stimulus
     */
    double GetIntracellularStimulus(double time);

    /**
     * Get the value of the intracellular stimulus.
     * This will always be in units of uA/cm^2.
     *
     * @param time  the time at which to evaluate the stimulus
     */
    double GetIntracellularAreaStimulus(double time);

    /**
     * Set whether this cell object exists in the context of a tissue simulation,
     * or can be used for single cell simulations.  This affects the units of the
     * intracellular stimulus (see GetIntracellularStimulus) and so is used by
     * GetIntracellularAreaStimulus to perform a units conversion if necessary.
     *
     * @param tissue  true if cell is in a tissue
     */
    void SetUsedInTissueSimulation(bool tissue=true);

    /**
     * Use CellML metadata to set up the default stimulus for this cell.
     * By default this method will always throw an exception.  For suitably annotated
     * models, PyCml will override this to provide a RegularStimulus as defined in
     * the CellML.
     */
    virtual void UseCellMLDefaultStimulus();

    /**
     *  [Ca_i] is needed for mechanics, so we explcitly have a Get method (rather than
     *  use a get by name type method, to avoid inefficiency when using different
     *  types of cells). This method by default throws an exception, so should be
     *  implemented in the concrete class if intracellular (cytosolic) calcium concentration is
     *  one of the state variables.
     */
    virtual double GetIntracellularCalciumConcentration();
    
    /**
     *  In electromechanics problems, the stretch is passed back to cell-model in case 
     *  mechano-electric feedback has been modelled. We define an empty method here.
     *  Stretch-dependent cell models should overload this method to use the input 
     *  stretch accordingly.
     *  @param stretch the stretch of the cell in the axial direction
     */
    virtual void SetStretch(double stretch)
    {
    }

    /**
     *  Empty method which can be over-ridden in concrete cell class which should
     *  go through the current state vector and go range checking on the values
     *  (eg check that concentrations are positive and gating variables are between
     *  zero and one). This method is called in the ComputeExceptVoltage method.
     */
    virtual void VerifyStateVariables()
    {
//// This code is for the future, but commented out at the moment due to the memory increase
//// it will introduce. See #794.
////
//// DOXYGEN DESCRIPTION NEEDS CHANGING ONCE THIS IS BROUGHT IN
////
////
//        for(std::set<unsigned>::iterator iter = mGatingVariableIndices.begin();
//            iter != mGatingVariableIndices.end();
//            ++iter)
//        {
//            double value = mStateVariables[*iter];
//            if(value<0.0)
//            {
//                std::stringstream error;
//                error << "State variable " << *iter << ", a gating variable, has gone negative";
//                EXCEPTION(DumpState(error.str()));
//            }
//            if(value>1.0)
//            {
//                std::stringstream error;
//                error << "State variable " << *iter << ", a gating variable, has become greater than one";
//                EXCEPTION(DumpState(error.str()));
//            }
//        }
//
//        for(std::set<unsigned>::iterator iter = mConcentrationIndices.begin();
//            iter != mConcentrationIndices.end();
//            ++iter)
//        {
//            if(mStateVariables[*iter] < 0.0)
//            {
//                std::stringstream error;
//                error << "State variable " << *iter << ", a concentration, has gone negative";
//                EXCEPTION(DumpState(error.str()));
//            }
//        }
    }



    ////////////////////////////////////////////////////////////////////////
    //  METHODS NEEDED BY FAST CARDIAC CELLS
    ////////////////////////////////////////////////////////////////////////

    /**
     * This should be implemented by fast/slow cardiac cell subclasses, and
     *  \li set the state
     *  \li initialise the cell
     *  \li \b SET #mNumberOfStateVariables \b CORRECTLY
     *      (as this would not have been known in the constructor.
     *
     * \note  This \e must be implemented by fast/slow cardiac cell subclasses.
     *
     * @param state  whether this cell is in fast or slow mode.
     */
    virtual void SetState(CellModelState state);

    /**
     * Set the slow variables. Should only be valid in fast mode.
     *
     * \note  This \e must be implemented by fast/slow cardiac cell subclasses.
     *
     * @param rSlowValues  values for the slow variables
     */
    virtual void SetSlowValues(const std::vector<double> &rSlowValues);

    /**
     * Get the current values of the slow variables. Should only be valid in slow mode.
     *
     * \note  This \e must be implemented by fast/slow cardiac cell subclasses.
     *
     * @param rSlowValues  will be filled in with the values of the slow variables on return.
     */
    virtual void GetSlowValues(std::vector<double>& rSlowValues);

    /** Get whether this cell is a fast or slow version.
     *
     * \note  This \e must be implemented by fast/slow cardiac cell subclasses.
     */
    virtual bool IsFastOnly();

    /**
     * In a multiscale simulation a cut-down cell model can be run:
     *  - fast values are calculated according to the CellML definition
     *  - slow values are interpolated on synchronisation time-steps.
     * There's a chance that linear interpolation/extrapolation may push
     * some gating variable out of the range [0, 1].  This method alters
     * any values which are out-of-range.
     *
     * \note  This \e must be implemented by fast/slow cardiac cell subclasses.
     *
     * @param rSlowValues A vector of the slow values for a particular cell after they have been interpolated from nearby coarse cells
     */
    virtual void AdjustOutOfRangeSlowValues(std::vector<double>& rSlowValues);

    /**
     * Get the number of slow variables for the cell model
     * (irrespective of whether in fast or slow mode).
     *
     * \note  This \e must be implemented by fast/slow cardiac cell subclasses.
     */
    virtual unsigned GetNumSlowValues();

    /**
     * For boost archiving use only
     * (that's why the 'consts' are required)
     *
     * @return The Intracellular stimulus pointer
     */
    const boost::shared_ptr<AbstractStimulusFunction> GetStimulusFunction() const;

    /**
     * For boost archiving use only
     * (that's why the 'consts' are required)
     *
     * @return pointer to the ODE solver being used
     */
    const boost::shared_ptr<AbstractIvpOdeSolver> GetSolver() const;

};

CLASS_IS_ABSTRACT(AbstractCardiacCell)
BOOST_CLASS_VERSION(AbstractCardiacCell, 1)

#endif /*ABSTRACTCARDIACCELL_HPP_*/
