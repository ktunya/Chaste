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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTUSINGCONTACTINHIBITIONCELLCYCLEMODELTUTORIAL_HPP_
#define TESTUSINGCONTACTINHIBITIONCELLCYCLEMODELTUTORIAL_HPP_

/*
 * = An example showing how to use the contact inhibition cell cycle model (with the contact inhibition simulator) =
 *
 * == Introduction ==
 *
 * In this tutorial, we will show how to use a simple implementation of the contact inhibition cell-cycle mode,
 * that stops cell division when the volume of the cell is smaller than a critical value.
 * We consider 2-D configurations for which cells are trapped in a square box. Firstly, all the cells are contact
 * inhibited and we study the effect of the critical volume upon the global cell density. Secondly, we consider
 * a mix of normal cells (contact inhibited) and tumour cells (not inhibited) and study the growth of the
 * tumour cells within the box.
 *
 *
 * == Including header files ==
 *
 * We begin by including the necessary header files. */
#include <cxxtest/TestSuite.h>

/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"
/* The next header include the NEVER_REACHED macro, used in one of the methods below. */
#include "Exception.hpp"

/*
 * The next header file defines the contact inhibition cell-cycle model that inherits from {{{AbstractSimpleCellCycleModel}}}.
 * The duration of the G1 phase depends on the deviation from a target volume (or area/length in 2D/1D): if the volume is
 * lower than a given fraction of the target volume, the G1 phase continues. The target volume and the critical fraction
 * are indicated in the user's Test file, and compared to the real volumes stored in {{{CellwiseData}}}, a singleton class.
 * This model allows for quiescence imposed by transient periods of high stress, followed by relaxation. Note that
 * in this cell cycle model, quiescence is implemented only by extending the G1 phase. Therefore, if a cell
 * is compressed during G2 or S phases then it will still divide, and thus cells whose volumes are smaller
 * than the given threshold may still divide.
 */
#include "ContactInhibitionCellCycleModel.hpp"

/*
 * The next header is the simulation class corresponding to the contact inhibition cell-cycle model. 
 * The essential difference with other simulators is that {{{CellWiseData}}} is updated with the 
 * volumes of the Voronoi elements representing each cell.
 */
#include "VolumeTrackedMeshBasedSimulation.hpp"
/* The remaining header files define classes that will be also be used and are presented in other tutorials. */
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SimulationTime.hpp"
#include "CellLabel.hpp"
#include "MutableMesh.hpp"
#include "PlaneBoundaryCondition.hpp"

/* We first define the global test class that inherits from {{{CxxTest::TestSuite}}}. */
class TestUsingContactInhibitionCellCycleModelTutorial : public CxxTest::TestSuite
{
public:
    /*
     * == Testing only normal cells ==
     *
     * We begin by testing the behaviour of normal cells (contact inhibited) trapped in a box.
     */
    void TestContactInhibitionInBox()
    {
        /* We must first set the start time. In addition, it is advisable to reset
         * the values of all model parameters. Recall that {{{SimulationTime}}} is
         * a ''singleton'' class; this means one and only one of each of this object
         * is instantiated at any time, and that single object is accessible from
         * anywhere in the code. As a result, we do not need to keep passing around
         * the present time. */
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        /* We use the honeycomb mesh generator to create a honeycomb mesh and
         * the associated mutable mesh. */
        HoneycombMeshGenerator generator(3, 3);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* We now create a vector of cell pointers. */
        std::vector<CellPtr> cells;

        /* We then define the mutation state of the cells we are working with. We will just consider
         * wild type mutations here. */
        MAKE_PTR(WildTypeCellMutationState, p_state);

        /* We now create a cell-cycle (only contact inhibited) model for these cells and loop over the
         * nodes of the mesh to create as many elements in the vector of cell pointers as there are
         * in the initial mesh. */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetCellProliferativeType(STEM);
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-(double)i);
            p_cycle_model->SetQuiescentVolumeFraction(0.5);
            p_cycle_model->SetEquilibriumVolume(1.0);
            p_cycle_model->SetStemCellG1Duration(0.5);
            p_cycle_model->SetTransitCellG1Duration(0.5);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        /* We now create a cell population, that takes several inputs: the mesh (for the position); and
         * the vector of cell pointers (for cycles and states)*/
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* Next, we create a force (springs) to be applied between cells and set up a
         * cut-off length beyond which cells stop interacting. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);

        /* To keep track of the volumes of the cells that are used in the contact inhibition cell-cycle,
         * we use the singleton class {{{CellWiseData}}}. Here, we just initialise it with one variable
         * and associate it with the cell population: */
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(1.0, cell_population.GetLocationIndexUsingCell(*cell_iter));
        }

        /*  Then, we define the contact {{{VolumeTrackedMeshBasedSimulation}}} class, that automatically updates the volumes of the cells
         * in {{{CellWiseData}}}. We also set up the output directory, the end time and pass the forces to the simulator.
         *
         */
        VolumeTrackedMeshBasedSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestBoxContactInhibition");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(50.0);
        simulator.AddForce(p_force);

        /*  To study the behaviour of the cells with varying volume, we trap them in a box, i.e., between
         *  4 plane boundary conditions. These planes are indicated by a point and a normal and then passed
         *  to the simulator:
         */
        c_vector<double,2> point = zero_vector<double>(2);
		c_vector<double,2> normal = zero_vector<double>(2);
		point(0) = 0.0;
		point(1) = 0.0;
		normal(0) = -1.0;
		normal(1) = 0.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // x>0
		simulator.AddCellPopulationBoundaryCondition(p_bc1);
		point(0) = 2.5;
		point(1) = 0.0;
		normal(0) = 1.0;
		normal(1) = 0.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // x<2.5
		simulator.AddCellPopulationBoundaryCondition(p_bc2);
		point(0) = 0.0;
		point(1) = 0.0;
		normal(0) = 0.0;
		normal(1) = -1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
		simulator.AddCellPopulationBoundaryCondition(p_bc3);
		point(0) = 0.0;
		point(1) = 2.5;
		normal(0) = 0.0;
		normal(1) = 1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<2.5
		simulator.AddCellPopulationBoundaryCondition(p_bc4);

        /* Test that the Solve() method does not throw any exceptions. */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* To make sure that the simulation is stopping with the correct cell volumes, we test that
         * the volume of the Voronoi elements associated with each cell is less than the critical
         * value for inhibition.
         */
//        cell_population.CreateVoronoiTessellation();
//        for (MeshBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
//             cell_iter != cell_population.End();
//             ++cell_iter)
//        {
//            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
//            TS_ASSERT_LESS_THAN(cell_population.GetVolumeOfVoronoiElement(node_index), 0.51);
//        }

        /* Finally, as in previous cell-based Chaste tutorials, we call {{{Destroy()}}} on the singleton classes. */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
    }

    /*
     * == Testing normal and tumour cells ==
     *
     * We now test the behaviour of a mixture of normal and tumour cells
     */
    void TestContactInhibitionInBoxWithMutants()
    {
        /* Set up SimulationTime. */
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        /* Create a simple mesh. */
        HoneycombMeshGenerator generator(2, 2, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* Create cell state. */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;

        /* Create cell cycle. The difference here is that one of the cell is not contact-inhibited, but rather
         * is defined by a stochastic duration generation-based cell-cycle model. */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            if (i==1)
            {
                StochasticDurationGenerationBasedCellCycleModel* p_cycle_model = new StochasticDurationGenerationBasedCellCycleModel();
                p_cycle_model->SetCellProliferativeType(TRANSIT);
                p_cycle_model->SetMaxTransitGenerations(UINT_MAX);
                p_cycle_model->SetTransitCellG1Duration(1);

                CellPtr p_cell(new Cell(p_state, p_cycle_model));
                p_cell->SetBirthTime(0.0);
                cells.push_back(p_cell);
            }
            else
            {
                ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
                p_cycle_model->SetCellProliferativeType(STEM);
                p_cycle_model->SetDimension(2);
                p_cycle_model->SetBirthTime(0.0);
                p_cycle_model->SetQuiescentVolumeFraction(0.8);
                p_cycle_model->SetEquilibriumVolume(1.0);
                p_cycle_model->SetStemCellG1Duration(0.5);

                CellPtr p_cell(new Cell(p_state, p_cycle_model));
                p_cell->InitialiseCellCycleModel();
                cells.push_back(p_cell);
            }
        }

        /* Create a cell population. */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* Create a force law. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);

        /* Create a singleton class to store the volume of the cells and initialise it. */
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        /* Create a contact inhibition simulator. */
        VolumeTrackedMeshBasedSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestBoxMixContactInhibition");
        simulator.SetEndTime(10.0);
        simulator.SetSamplingTimestepMultiple(12);
        simulator.AddForce(p_force);

        /* Create some boundary conditions and pass them to the simulation. */
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        point(0) = -4.0;
        point(1) = 0.0;
        normal(0) = -1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        point(0) = 4.0;
        point(1) = 0.0;
        normal(0) = 1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<2
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
        point(0) = 0.0;
        point(1) = -4.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y<2
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        point(0) = 0.0;
        point(1) = 4.0;
        normal(0) = 0.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<2
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        /* Test that the Solve() method does not throw any exceptions. */
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        /* Finally, as in previous cell-based Chaste tutorials, we call {{{Destroy()}}} on the singleton classes. */
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
    }
};
#endif /*TESTUSINGCONTACTINHIBITIONCELLCYCLEMODELTUTORIAL_HPP_*/
