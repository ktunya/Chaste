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

#include "OffLatticeSimulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"

template<unsigned DIM>
OffLatticeSimulation<DIM>::OffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                bool deleteCellPopulationInDestructor,
                                                bool initialiseCells)
    : AbstractCellBasedSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
      mOutputNodeVelocities(false)
{
    if (!dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OnLatticeSimulations require a subclass of AbstractOffLatticeCellPopulation.");
    }

    // Different time steps are used for cell-centre and vertex-based simulations
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        this->mDt = 1.0/120.0; // 30 seconds
    }
    else
    {
        this->mDt = 0.002; // smaller time step required for convergence/stability
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::AddForce(boost::shared_ptr<AbstractForce<DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<DIM> > pBoundaryCondition)
{
    mBoundaryConditions.push_back(pBoundaryCondition);
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::UpdateCellLocationsAndTopology()
{
    // Calculate forces
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);

    // Initialise a vector of forces on node
    std::vector<c_vector<double, DIM> > forces(this->mrCellPopulation.GetNumNodes(), zero_vector<double>(DIM));

/**\todo Is it faster to preallocate and have forces as a member variable? see #1890**/
//    // First set all the forces to zero
//    for (unsigned i=0; i<forces.size(); i++)
//    {
//         forces[i].clear();
//    }
//
//    // Then resize the std::vector if the number of cells has increased or decreased
//    // (note this should be done after the above zeroing)
//    unsigned num_nodes = mrCellPopulation.GetNumNodes();
//    if (num_nodes != forces.size())
//    {
//        forces.resize(num_nodes, zero_vector<double>(DIM));
//    }

    // Now add force contributions from each AbstractForce
    for (typename std::vector<boost::shared_ptr<AbstractForce<DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        (*iter)->AddForceContribution(forces, this->mrCellPopulation);
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

    // Update node positions
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
    UpdateNodePositions(forces);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
}

template<unsigned DIM>
c_vector<double, DIM> OffLatticeSimulation<DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    /**
     * \todo Could remove this dynamic_cast by moving the code block below into
     * AbstractCentreBasedCellPopulation::AddCell(), allowing it to be overruled by
     * this method when overridden in subclasses. See also comment on #1093.
     */
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        // Location of parent and daughter cells
        c_vector<double, DIM> parent_coords = this->mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, DIM> daughter_coords;

        // Get separation parameter
        double separation = static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->GetMeinekeDivisionSeparation();

        // Make a random direction vector of the required length
        c_vector<double, DIM> random_vector;

        /*
         * Pick a random direction and move the parent cell backwards by 0.5*separation
         * in that direction and return the position of the daughter cell 0.5*separation
         * forwards in that direction.
         */
        switch (DIM)
        {
            case 1:
            {
                double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

                random_vector(0) = 0.5*separation*random_direction;
                break;
            }
            case 2:
            {
                double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

                random_vector(0) = 0.5*separation*cos(random_angle);
                random_vector(1) = 0.5*separation*sin(random_angle);
                break;
            }
            case 3:
            {
                double random_zenith_angle = M_PI*RandomNumberGenerator::Instance()->ranf(); // phi
                double random_azimuth_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf(); // theta

                random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(2) = 0.5*separation*cos(random_zenith_angle);
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }

        parent_coords = parent_coords - random_vector;
        daughter_coords = parent_coords + random_vector;

        // Set the parent to use this location
        ChastePoint<DIM> parent_coords_point(parent_coords);
        unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        this->mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else
    {
        ///\todo do something for vertex models here
        return zero_vector<double>(DIM);
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::WriteVisualizerSetupFile()
{
    for (unsigned i=0; i<this->mForceCollection.size(); i++)
    {
        boost::shared_ptr<AbstractForce<DIM> > p_force = this->mForceCollection[i];
        if (boost::dynamic_pointer_cast<AbstractTwoBodyInteractionForce<DIM> >(p_force))
        {
            double cutoff = (boost::static_pointer_cast<AbstractTwoBodyInteractionForce<DIM> >(p_force))->GetCutOffLength();
            *(this->mpVizSetupFile) << "Cutoff\t" << cutoff << "\n";
        }
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rNodeForces)
{
    unsigned num_nodes = this->mrCellPopulation.GetNumNodes();

    /*
     * Get the previous node positions (these may be needed when applying boundary conditions,
     * e.g. in the case of immotile cells)
     */
    std::vector<c_vector<double, DIM> > old_node_locations;
    old_node_locations.reserve(num_nodes);
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        old_node_locations[node_index] = this->mrCellPopulation.GetNode(node_index)->rGetLocation();
    }

    // Update node locations
    static_cast<AbstractOffLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->UpdateNodeLocations(rNodeForces, this->mDt);

    // Apply any boundary conditions
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        (*bcs_iter)->ImposeBoundaryCondition(old_node_locations);
    }

    // Verify that each boundary condition is now satisfied
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        if (!((*bcs_iter)->VerifyBoundaryCondition()))
        {
            EXCEPTION("The cell population boundary conditions are incompatible.");
        }
    }

    // Write node velocities to file if required
    if (mOutputNodeVelocities)
    {
        if (SimulationTime::Instance()->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0)
        {
            *mpNodeVelocitiesFile << SimulationTime::Instance()->GetTime() << "\t";
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                // We should never encounter deleted nodes due to where this method is called by Solve()
                assert(!this->mrCellPopulation.GetNode(node_index)->IsDeleted());

                // Check that results should be written for this node
                bool is_real_node = true;

                if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&this->mrCellPopulation))
                {
                    if (static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->IsGhostNode(node_index))
                    {
                        // If this node is a ghost node then don't record its velocity
                        is_real_node = false;
                    }
                    else
                    {
                        // We should never encounter nodes associated with dead cells due to where this method is called by Solve()
                        assert(!this->mrCellPopulation.GetCellUsingLocationIndex(node_index)->IsDead());
                    }
                }

                // Write node data to file
                if (is_real_node)
                {
                    const c_vector<double,DIM>& position = this->mrCellPopulation.GetNode(node_index)->rGetLocation();
                    double damping_constant = static_cast<AbstractOffLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->GetDampingConstant(node_index);
                    c_vector<double, DIM> velocity = this->mDt * rNodeForces[node_index] / damping_constant;

                    *mpNodeVelocitiesFile << node_index  << " ";
                    for (unsigned i=0; i<DIM; i++)
                    {
                        *mpNodeVelocitiesFile << position[i] << " ";
                    }
                    for (unsigned i=0; i<DIM; i++)
                    {
                        *mpNodeVelocitiesFile << velocity[i] << " ";
                    }
                }
            }
            *mpNodeVelocitiesFile << "\n";
        }
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::SetupSolve()
{
    // First call method on base class
    AbstractCellBasedSimulation<DIM>::SetupSolve();

    if (mOutputNodeVelocities)
    {
        OutputFileHandler output_file_handler2(this->mSimulationOutputDirectory+"/", false);
        mpNodeVelocitiesFile = output_file_handler2.OpenOutputFile("nodevelocities.dat");
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::AfterSolve()
{
    // First call method on base class
    AbstractCellBasedSimulation<DIM>::AfterSolve();

    if (mOutputNodeVelocities)
    {
        mpNodeVelocitiesFile->close();
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over forces
    *rParamsFile << "\n\t<Forces>\n";
    for (typename std::vector<boost::shared_ptr<AbstractForce<DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        // Output force details
        (*iter)->OutputForceInfo(rParamsFile);
    }
    *rParamsFile << "\t</Forces>\n";

    // Loop over cell population boundary conditions
    *rParamsFile << "\n\t<CellPopulationBoundaryConditions>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<DIM> > >::iterator iter = mBoundaryConditions.begin();
         iter != mBoundaryConditions.end();
         ++iter)
    {
        // Output cell Boundary condition details
        (*iter)->OutputCellPopulationBoundaryConditionInfo(rParamsFile);
    }
    *rParamsFile << "\t</CellPopulationBoundaryConditions>\n";
}

template<unsigned DIM>
bool OffLatticeSimulation<DIM>::GetOutputNodeVelocities()
{
    return mOutputNodeVelocities;
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::SetOutputNodeVelocities(bool outputNodeVelocities)
{
    mOutputNodeVelocities = outputNodeVelocities;
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<OutputNodeVelocities>" << mOutputNodeVelocities << "</OutputNodeVelocities>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OffLatticeSimulation<1>;
template class OffLatticeSimulation<2>;
template class OffLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeSimulation)