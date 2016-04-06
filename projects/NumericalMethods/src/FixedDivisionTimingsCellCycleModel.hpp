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

#ifndef FIXEDDIVISIONTIMINGSCELLCYCLEMODEL_HPP_
#define FIXEDDIVISIONTIMINGSCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellCycleModel.hpp"

/**
 * This class will ensure that all cells divide simultaneously at a precisely specified time.
 *
 * Useful for studying the differences between several simulations while controlling for division 
 * timings.
 */
class FixedDivisionTimingsCellCycleModel : public AbstractCellCycleModel
{

private:

	/** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
    }

public:

    double targetTimeBetweenDivisions; 

    FixedDivisionTimingsCellCycleModel(double timeBetweenDivisions);


    virtual ~FixedDivisionTimingsCellCycleModel();


    virtual void UpdateCellCyclePhase();


    virtual AbstractCellCycleModel* CreateCellCycleModel();


    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(FixedDivisionTimingsCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a FixedDivisionTimingsCellCycleModel.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const FixedDivisionTimingsCellCycleModel* t, const unsigned int file_version)
{
  ar << t->targetTimeBetweenDivisions;
}

/**
 * De-serialize constructor parameters and initialize a FixedDivisionTimingsCellCycleModel.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, FixedDivisionTimingsCellCycleModel* t, const unsigned int file_version)
{
    // Retrieve other member variables
    double targetDivisionInterval;
    ar >> targetDivisionInterval;

    // Invoke inplace constructor to initialise instance
    ::new(t)FixedDivisionTimingsCellCycleModel(targetDivisionInterval);
}
}
} // namespace ...


#endif /*FIXEDDIVISIONTIMINGSCELLCYCLEMODEL_HPP_*/