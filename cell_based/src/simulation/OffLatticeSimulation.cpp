/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "OffLatticeSimulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "T2SwapCellKiller.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"
#include "Cylindrical2dMesh.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "StepSizeException.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OffLatticeSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                bool deleteCellPopulationInDestructor,
                                                bool initialiseCells,
                                                boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > numericalMethod,
                                                bool isAdaptive)
    : AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
      mNumericalMethod(numericalMethod),
      mAdaptive(isAdaptive)
{

    if (!dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OffLatticeSimulations require a subclass of AbstractOffLatticeCellPopulation.");
    }

    // Different time steps are used for cell-centre and vertex-based simulations
    if (bool(dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    {
        this->mDt = 1.0/120.0; // 30 seconds
    }
    else if (bool(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        this->mDt = 0.002; // smaller time step required for convergence/stability

        // For VertexBasedCellPopulations we automatically add a T2SwapCellKiller. In order to inhibit T2 swaps
        // the user needs to set the threshold for T2 swaps in the mesh to 0.
        VertexBasedCellPopulation<SPACE_DIM>* p_vertex_based_cell_population = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
        MAKE_PTR_ARGS(T2SwapCellKiller<SPACE_DIM>, T2_swap_cell_killer, (p_vertex_based_cell_population));
        this->AddCellKiller(T2_swap_cell_killer);
    }
    else
    {
        // All classes derived from AbstractOffLatticeCellPopulation are covered by the above (except user-derived classes),
        // i.e. if you want to use this method with your own subclass of AbstractOffLatticeCellPopulation, then simply
        // comment out the line below
        NEVER_REACHED;
    }

    mNumericalMethod->SetCellPopulation(dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation));
    mNumericalMethod->SetForceCollection(&mForceCollection);
    mNumericalMethod->SetAdaptive(&mAdaptive);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::AddForce(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RemoveAllForces()
{
    mForceCollection.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > pBoundaryCondition)
{
    mBoundaryConditions.push_back(pBoundaryCondition);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RemoveAllCellPopulationBoundaryConditions()
{
    mBoundaryConditions.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::GetNumericalMethod() const{
    return mNumericalMethod;
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const bool OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::IsAdaptive() const{
    return mAdaptive;
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::SetAdaptive(bool isAdaptive){
    mAdaptive = isAdaptive;
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateCellLocationsAndTopology()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);

    double timeAdvancedSoFar = 0; 
    double targetTimeStep  = this->mDt;
    double currentTimeStep = this->mDt;

    while(timeAdvancedSoFar < targetTimeStep){

        // Store the initial node positions (these may be needed when applying boundary conditions)    
        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations;

        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
            node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
        {
            old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

        // Try to update node positions according to the numerical method 
        try{

            mNumericalMethod->UpdateAllNodePositions(currentTimeStep);
            ApplyBoundaries(old_node_locations);

            // Successful timestep! Update timeAdvancedSoFar and increase the currentTimeStep
            // (by 1% for now)
            timeAdvancedSoFar += currentTimeStep;
            if(mAdaptive){
                currentTimeStep = fmin(1.01*currentTimeStep, targetTimeStep - timeAdvancedSoFar);
            }

        }catch(StepSizeException* e){
            // Detects if a node has travelled too far in a single time step
            if(mAdaptive){
                // If adaptivity is switched on, revert node locations and choose a suitable
                // smaller time step
                RevertToOldLocations(old_node_locations);
                currentTimeStep = fmin(e->suggestedNewStep, targetTimeStep - timeAdvancedSoFar); 
            }else{
                // If adaptivity is switched off, terminate with an error
                EXCEPTION(e->what());
            }
        }

    }

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RevertToOldLocations(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations){
    
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
        ++node_iter)
    {
        (node_iter)->rGetModifiableLocation() = old_node_locations[&(*node_iter)];
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::ApplyBoundaries(std::map<Node<SPACE_DIM>*,c_vector<double, SPACE_DIM> > old_node_locations){

    // Apply any boundary conditions
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        (*bcs_iter)->ImposeBoundaryCondition(old_node_locations);
    }

    // Verify that each boundary condition is now satisfied
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        if (!((*bcs_iter)->VerifyBoundaryCondition()))
        {
            EXCEPTION("The cell population boundary conditions are incompatible.");
        }
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    /**
     * \todo #2400 Could remove this dynamic_cast by moving the code block below into
     * AbstractCentreBasedCellPopulation::AddCell(), allowing it to be overruled by
     * this method when overridden in subclasses. See also comment on #1093.
     */
    // If it is not vertex based...
    if (bool(dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))))
    {
        // Location of parent and daughter cells
        c_vector<double, SPACE_DIM> parent_coords = this->mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, SPACE_DIM> daughter_coords;

        // Get separation parameter
        double separation = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))->GetMeinekeDivisionSeparation();

        // Make a random direction vector of the required length
        c_vector<double, SPACE_DIM> random_vector;

        /*
         * Pick a random direction and move the parent cell backwards by 0.5*separation
         * in that direction and return the position of the daughter cell 0.5*separation
         * forwards in that direction.
         */
        switch (SPACE_DIM)
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
                /*
                 * Note that to pick a random point on the surface of a sphere, it is incorrect
                 * to select spherical coordinates from uniform distributions on [0, 2*pi) and
                 * [0, pi) respectively, since points picked in this way will be 'bunched' near
                 * the poles. See #2230.
                 */
                double u = RandomNumberGenerator::Instance()->ranf();
                double v = RandomNumberGenerator::Instance()->ranf();

                double random_azimuth_angle = 2*M_PI*u;
                double random_zenith_angle = std::acos(2*v - 1);

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
        ChastePoint<SPACE_DIM> parent_coords_point(parent_coords);
        unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        this->mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else if (bool(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation))))
    {
        VertexBasedCellPopulation<SPACE_DIM>* p_vertex_population = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation));
        boost::shared_ptr<AbstractVertexBasedDivisionRule<SPACE_DIM> > p_division_rule = p_vertex_population->GetVertexBasedDivisionRule();

        return p_division_rule->CalculateCellDivisionVector(pParentCell, *p_vertex_population);
    }
    else
    {
        // All classes derived from AbstractOffLatticeCellPopulation are covered by the above (except user-derived classes),
        // i.e. if you want to use this class with your own subclass of AbstractOffLatticeCellPopulation, then simply
        // comment out the line below
        NEVER_REACHED;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::WriteVisualizerSetupFile()
{
    if (PetscTools::AmMaster())
    {
        for (unsigned i=0; i<this->mForceCollection.size(); i++)
        {
            // This may cause compilation problems, probably due to AbstractTwoBodyInteractionForce not having two template parameters
            ///\todo Check whether this comment is still valid

            boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > p_force = this->mForceCollection[i];
            if (boost::dynamic_pointer_cast<AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM> >(p_force))
            {
                double cutoff = (boost::static_pointer_cast<AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM> >(p_force))->GetCutOffLength();
                *(this->mpVizSetupFile) << "Cutoff\t" << cutoff << "\n";
            }
        }

        // This is a quick and dirty check to see if the mesh is periodic
        if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&this->mrCellPopulation)))
        {
           if (bool(dynamic_cast<Cylindrical2dMesh*>(&(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))->rGetMesh()))))
           {
               *this->mpVizSetupFile << "MeshWidth\t" << this->mrCellPopulation.GetWidth(0) << "\n";
           }
        }
        else if (bool(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&this->mrCellPopulation)))
        {
           if (bool(dynamic_cast<Cylindrical2dVertexMesh*>(&(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation))->rGetMesh()))))
           {
               *this->mpVizSetupFile << "MeshWidth\t" << this->mrCellPopulation.GetWidth(0) << "\n";
           }
        }
    }
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::SetupSolve()
{
    // Clear all forces
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over forces
    *rParamsFile << "\n\t<Forces>\n";
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        // Output force details
        (*iter)->OutputForceInfo(rParamsFile);
    }
    *rParamsFile << "\t</Forces>\n";

    // Loop over cell population boundary conditions
    *rParamsFile << "\n\t<CellPopulationBoundaryConditions>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mBoundaryConditions.begin();
         iter != mBoundaryConditions.end();
         ++iter)
    {
        // Output cell Boundary condition details
        (*iter)->OutputCellPopulationBoundaryConditionInfo(rParamsFile);
    }
    *rParamsFile << "\t</CellPopulationBoundaryConditions>\n";

    // Output numerical method details
    *rParamsFile << "\n\t<NumericalMethod>\n";
    mNumericalMethod->OutputNumericalMethodInfo(rParamsFile);
    *rParamsFile << "\t\t<AdaptiveStepSize>" << (int)mAdaptive << "</AdaptiveStepSize>\n";
    *rParamsFile << "\t</NumericalMethod>\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(rParamsFile);
}

///////// Explicit instantiation
template class OffLatticeSimulation<1,1>;
template class OffLatticeSimulation<1,2>;
template class OffLatticeSimulation<2,2>;
template class OffLatticeSimulation<1,3>;
template class OffLatticeSimulation<2,3>;
template class OffLatticeSimulation<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulation)
