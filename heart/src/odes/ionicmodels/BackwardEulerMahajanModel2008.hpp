#ifndef BackwardEulerMahajanModel2008_HPP_
#define BackwardEulerMahajanModel2008_HPP_

//! @file
//!
//! This source file was generated from CellML.
//!
//! Model: mahajan_shiferaw_2008_version01
//!
//! Processed by pycml - CellML Tools in Python
//!     (translate: 8196, pycml: 8196)
//! on Wed Mar 10 17:46:04 2010
//!
//! <autogenerated>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractBackwardEulerCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"

class BackwardEulerMahajanModel2008 : public AbstractBackwardEulerCardiacCell<15>
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractBackwardEulerCardiacCell<15> >(*this);
    }


public:
    BackwardEulerMahajanModel2008(boost::shared_ptr<AbstractIvpOdeSolver> /* unused; should be empty */, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~BackwardEulerMahajanModel2008();
    void VerifyStateVariables();
    double GetIIonic();
    void ComputeResidual(const double rCurrentGuess[15], double rResidual[15]);
    void ComputeJacobian(const double rCurrentGuess[15], double rJacobian[15][15]);
    void UpdateTransmembranePotential(double var_Environment__time);
    void ComputeOneStepExceptVoltage(double var_Environment__time);
};


// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BackwardEulerMahajanModel2008)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const BackwardEulerMahajanModel2008 * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, BackwardEulerMahajanModel2008 * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)BackwardEulerMahajanModel2008(p_solver, p_stimulus);
        }

    }

}

#endif // BackwardEulerMahajanModel2008_HPP_
