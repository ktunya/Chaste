/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef CRYPTSIMULATION2DBOUNDARYCONDITION_HPP_
#define CRYPTSIMULATION2DBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A boundary condition class for use with CryptSimulation2d, which pins stem cells
 * in the absence of a Wnt stimulus, and optionally prevents cells moving below the
 * y=0 boundary via random jiggling.
 */
class CryptSimulation2dBoundaryCondition : public AbstractCellPopulationBoundaryCondition<2>
{
private:

    /**
     * Whether to jiggle the cells on the bottom surface, initialised to false
     * in the constructor.
     */
    bool mUseJiggledBottomCells;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<2> >(*this);
        archive & mUseJiggledBottomCells;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     */
    CryptSimulation2dBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation);

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::vector< c_vector<double, 2> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Set method for mUseJiggledBottomCells
     *
     * @param useJiggledBottomCells whether to jiggle the cells on the bottom surface
     */
    void SetUseJiggledBottomCells(bool useJiggledBottomCells);

    /** Get method for mUseJiggledBottomCells. */
    bool GetUseJiggledBottomCells();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2dBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation2dBoundaryCondition.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation2dBoundaryCondition * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialize a CryptSimulation2dBoundaryCondition.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation2dBoundaryCondition * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSimulation2dBoundaryCondition(p_cell_population);
}
}
} // namespace ...

#endif /* CRYPTSIMULATION2DBOUNDARYCONDITION_HPP_ */
