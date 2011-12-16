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

#ifndef TESTCRYPTSIMULATION2DWITHMESHBASEDCELLPOPULATION_HPP_
#define TESTCRYPTSIMULATION2DWITHMESHBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "TargetedCellKiller.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "SimpleWntCellCycleModel.hpp"

class TestCryptSimulation2dWithMeshBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    /**
     * Compare two meshes to see if they are 'the same'.  Doesn't check everything,
     * but is fairly thorough.  Used for testing serialization.
     */
    template<unsigned DIM>
    void CompareMeshes(MutableMesh<DIM,DIM>* pMesh1,
                       MutableMesh<DIM,DIM>* pMesh2)
    {
        TS_ASSERT_EQUALS(pMesh1->GetNumAllNodes(), pMesh2->GetNumAllNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumNodes(), pMesh2->GetNumNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumBoundaryNodes(), pMesh2->GetNumBoundaryNodes());

        for (unsigned i=0; i<pMesh1->GetNumAllNodes(); i++)
        {
            Node<DIM>* p_node = pMesh1->GetNode(i);
            Node<DIM>* p_node2 = pMesh2->GetNode(i);
            TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
            TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
            TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
            for (unsigned j=0; j<DIM; j++)
            {
                TS_ASSERT_DELTA(p_node->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-16);
            }
        }

        TS_ASSERT_EQUALS(pMesh1->GetNumElements(), pMesh2->GetNumElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumAllElements(), pMesh2->GetNumAllElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumBoundaryElements(), pMesh2->GetNumBoundaryElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumAllBoundaryElements(), pMesh2->GetNumAllBoundaryElements());

        typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator iter2 = pMesh2->GetElementIteratorBegin();
        for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator iter = pMesh1->GetElementIteratorBegin();
             iter != pMesh1->GetElementIteratorEnd();
             ++iter, ++iter2)
        {
            TS_ASSERT_EQUALS(iter->GetNumNodes(), iter2->GetNumNodes());
            for (unsigned i=0; i<iter->GetNumNodes(); i++)
            {
                TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(i), iter2->GetNodeGlobalIndex(i));
            }
        }
    }

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestCryptSimulation2dExceptions()
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Try to set up off lattice simulation
        TS_ASSERT_THROWS_THIS(CryptSimulation2d simulator(node_based_cell_population),
            "CryptSimulation2d is to be used with MeshBasedCellPopulation or VertexBasedCellPopulation (or subclasses) only");
    }

    /**
     * Test the spring system.
     *
     * The cells in this test are given an intial age of 2.0 so that their
     * springs are at their natural length (i.e. we set birth time=-2.0).
     *
     * The mesh is initially a set of 10 by 10 squares, each square made up
     * of two triangles. The horizontal and vertical edges (springs) are at
     * rest length, the diagonals are longer, so this means the mesh skews
     * to a (sloughed) parallelogram, with each triangle trying to become
     * equilateral.
     *
     * If you want to view the results visually, set the end time to 24.0,
     * and the spring system will resemble a parallelogram. However we keep
     * the simulation time at 1.0 in order to keep the test short.
     */
    void Test2DSpringSystem() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;

        cells_generator.Generate(cells, &mesh, std::vector<unsigned>(), false, 0.0, 3.0, 6.5, 8.0);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        simulator.SetEndTime(1.0);

        TS_ASSERT_THROWS_THIS(simulator.Solve(),"OutputDirectory not set");

        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        simulator.SetOutputDirectory("Crypt2DSprings");
        simulator.SetEndTime(1.0);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetNoBirth(true);

        // Destroy the simulation time class because of failed solve
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        simulator.Solve();
        std::vector<double> node_0_location = simulator.GetNodeLocation(0);
        TS_ASSERT_DELTA(node_0_location[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(node_0_location[1], 0.0, 1e-12);

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DSprings", false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT_EQUALS(system(("diff " + node_results_file + " crypt/test/data/Crypt2DSpringsResults/results.viznodes").c_str()), 0);

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        TS_ASSERT_EQUALS(system(("diff " + cell_type_results_file + " crypt/test/data/Crypt2DSpringsResults/results.vizcelltypes").c_str()), 0);

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        TS_ASSERT_EQUALS(system(("diff " + elem_results_file + " crypt/test/data/Crypt2DSpringsResults/results.vizelements").c_str()), 0);
    }

    void TestWithMultipleCellKillers() throw (Exception)
    {
       unsigned cells_across = 7;
       unsigned cells_up = 11;
       unsigned thickness_of_ghost_layer = 4;

       HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
       MutableMesh<2,2>* p_mesh = generator.GetMesh();
       std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

       // Set up cells
       std::vector<CellPtr> cells;
       CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
       cells_generator.Generate(cells, p_mesh, location_indices, true);

       // Create cell population
       MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

       // Create crypt simulation from cell population and force law
       CryptSimulation2d simulator(crypt);
       simulator.SetOutputDirectory("CryptWithMultipleCellKillers");

       // Create a force law and pass it to the simulation
       MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
       simulator.AddForce(p_linear_force);

       // Create cell killer and pass in to crypt simulation.
       // They kill the first and second available cell,
       // which are attached to nodes 64 and 65 respectively.
       MAKE_PTR_ARGS(TargetedCellKiller<2>, p_killer1, (&crypt, 64));
       MAKE_PTR_ARGS(TargetedCellKiller<2>, p_killer2, (&crypt, 65));
       simulator.AddCellKiller(p_killer1);
       simulator.AddCellKiller(p_killer2);

       unsigned num_cells = crypt.GetNumRealCells();

       std::vector<bool> ghost_node_indices_before = crypt.rGetGhostNodes();
       unsigned num_ghosts_before = 0;
       for (unsigned i=0; i<ghost_node_indices_before.size(); i++)
       {
           if (ghost_node_indices_before[i])
           {
               num_ghosts_before++;
           }
       }

       // Just enough time to kill off the cells (and watch ghost mesh use force laws)
       double dt = 0.01;

       simulator.SetDt(dt);
       simulator.SetEndTime(dt);

       // Run simulation
       simulator.Solve();

       std::vector<bool> ghost_node_indices_after = crypt.rGetGhostNodes();
       unsigned num_ghosts_after = 0;
       for (unsigned i=0; i<ghost_node_indices_after.size(); i++)
       {
           if (ghost_node_indices_after[i])
           {
               num_ghosts_after++;
           }
       }

       // Check no nodes have been created
       TS_ASSERT_EQUALS(num_ghosts_after, num_ghosts_after);

       // All cells should have been removed in this time
       TS_ASSERT_EQUALS(crypt.GetNumRealCells(), num_cells-2u);
    }

    void TestUpdatePositions() throw (Exception)
    {
        // Create mesh
        HoneycombMeshGenerator generator(3, 3, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(cell_population);

        // No force law passed to simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        std::vector<c_vector<double, 2> > old_posns(p_mesh->GetNumNodes());
        std::vector<c_vector<double, 2> > forces(p_mesh->GetNumNodes());

        // Make up some forces
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            old_posns[i][0] = p_mesh->GetNode(i)->rGetLocation()[0];
            old_posns[i][1] = p_mesh->GetNode(i)->rGetLocation()[1];

            forces[i][0] = i*0.01;
            forces[i][1] = 2*i*0.01;
        }

        simulator.SetDt(0.01);
        simulator.UpdateNodePositions(forces);

        // Create a set of node indices corresponding to ghost nodes
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices_set;
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            node_indices.insert(p_mesh->GetNode(i)->GetIndex());
        }
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            location_indices_set.insert(location_indices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices_set.begin(), location_indices_set.end(),
                            std::inserter(ghost_node_indices, ghost_node_indices.begin()));

        for (AbstractCellPopulation<2>::Iterator cell_iter=simulator.rGetCellPopulation().Begin();
             cell_iter!=simulator.rGetCellPopulation().End();
             ++cell_iter)
        {
            unsigned index = simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 2> cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);

            if (old_posns[index][1] == 0) // stem
            {
                // No Wnt so shouldn't have been moved
                TS_ASSERT_DELTA(cell_location[0], old_posns[index][0], 1e-9);
                TS_ASSERT_DELTA(cell_location[1], old_posns[index][1], 1e-9);
            }
            else
            {
                TS_ASSERT_DELTA(cell_location[0], old_posns[index][0] +   index*0.01*0.01, 1e-9);
                TS_ASSERT_DELTA(cell_location[1], old_posns[index][1] + 2*index*0.01*0.01, 1e-9);
            }
        }
    }

    void Test2DCylindrical() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        double crypt_length = DBL_MAX; // so no sloughing
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);// true = mature cells

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.SetOutputCellMutationStates(true);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt);
        simulator.SetEndTime(0.1);

        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetUpdateCellPopulationRule(true);
        simulator.SetNoBirth(false);
        simulator.SetOutputDirectory("Crypt2DCylindrical");

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        std::vector<unsigned> cell_mutation_states = crypt.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states[0], crypt.GetNumRealCells());

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_nodes, ghost_cells.size());
        TS_ASSERT_EQUALS(number_of_cells, cells_across*cells_up+1u);  // 6 cells in a row*12 rows + 1 birth
        TS_ASSERT_EQUALS(number_of_nodes, number_of_cells+thickness_of_ghost_layer*2*cells_across);

        // Coverage of exceptions (after main test to avoid problems with SimulationTime).
        simulator.SetEndTime(10.0);
        simulator.SetOutputDirectory("");
        TS_ASSERT_THROWS_THIS(simulator.Solve(),"OutputDirectory not set");
        CellBasedEventHandler::Reset(); // otherwise event handler left in bad state after throw
    }


    void Test2DCylindricalMultipleDivisions() throw (Exception)
    {
        // Create a log of this test
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(2, "Crypt2DCylindricalMultipleDivisions");

        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 8;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        double crypt_length = cells_up*(sqrt(3)/2)*crypt_width/cells_across;

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true); // true = mature cells

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-11.5);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DCylindricalMultipleDivisions");
        simulator.SetEndTime(0.6);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Find the height of the current crypt
        double height_after_division = p_mesh->GetWidth(1);

        // Reset end time and run simulation
        simulator.SetEndTime(0.8);
        simulator.Solve();

        // Find the height of the current crypt
        double height_after_relaxation = p_mesh->GetWidth(1);

        TS_ASSERT_LESS_THAN(height_after_division, height_after_relaxation);

        // Reset end time and run simulation
        simulator.SetEndTime(2.0);
        simulator.Solve();

        // All fully differentiated cells have sloughed off
        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
             TS_ASSERT(cell_iter->GetCellCycleModel()->GetCellProliferativeType() != DIFFERENTIATED);
        }

        // Close the log file opened in this test
        LogFile::Close();
    }

    /*
     * This test compares the visualizer output from the previous test, Test2DCylindricalMultipleDivisions,
     * with a known file.
     *
     * The results of this should be a yellow crypt with a line of
     * blue stem cells at the base, and four labelled cells of different
     * colours in the centre. A couple of cells divide, the crypt stays
     * periodic and a couple of cells swap sides.
     *
     * Note - if the previous test is changed we need to update the file this test refers to.
     */
    void TestVisualizerOutput() throw (Exception)
    {
        // Work out where one of the previous tests wrote its files
        OutputFileHandler handler("Crypt2DCylindricalMultipleDivisions", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_ele(results_dir + "/results.vizelements", "crypt/test/data/Crypt2DCylindricalMultipleDivisions/results.vizelements");
        TS_ASSERT(comp_ele.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizelements crypt/test/data/Crypt2DCylindricalMultipleDivisions/results.vizelements").c_str()), 0);

        NumericFileComparison comp_nodes(results_dir + "/results.viznodes", "crypt/test/data/Crypt2DCylindricalMultipleDivisions/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles(1e-15));

        NumericFileComparison comp_celltypes(results_dir + "/results.vizcelltypes", "crypt/test/data/Crypt2DCylindricalMultipleDivisions/results.vizcelltypes");
        TS_ASSERT(comp_celltypes.CompareFiles(1e-15));

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizsetup crypt/test/data/Crypt2DCylindricalMultipleDivisions/results.vizsetup").c_str()), 0);
    }


    // This is a rubbish test - all cells start at birthTime = 0.
    // So bizarrely the crypt shrinks as the rest lengths are shortened!
    // But at least it uses Wnt cell cycle and runs reasonably quickly...
    // For a better test with more randomly distributed cell ages see the Nightly test pack.
    void TestWithWntDependentCells() throw (Exception)
    {
        // Create mesh
        double crypt_length = 22.0;
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), false);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicWnt");
        simulator.SetEndTime(0.3);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        std::vector<double> node_35_location = simulator.GetNodeLocation(35);

        TS_ASSERT_DELTA(node_35_location[0], 5.5000, 1e-4);

        // Old version of this test had cells with age zero, therefore small spring lengths.
        // Variable spring lengths now only associated with cell division.
        TS_ASSERT_DELTA(node_35_location[1], 4.33013, 1e-4);

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    // A better check that the loaded mesh is the same as that saved
    void TestMeshSurvivesSaveLoad() throw (Exception)
    {
        // Create mesh
        double crypt_length = 22.0;
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, false);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DMeshArchive");
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Memory leak (unconditional jump) without the following line.
        // The archiver assumes that a Solve has been called and simulation time has been set up properly.
        // In this test it hasn't so we need this to avoid memory leak.
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.1, 100);

        // And record current state of mesh
        MutableMesh<2,2>& r_mesh = (static_cast<MeshBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->rGetMesh();

        // Save
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        // Load
        CryptSimulation2d* p_simulator;
        p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("Crypt2DMeshArchive", 0);

        // Get the loaded mesh
        MutableMesh<2,2>& r_mesh2 = (static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->rGetMesh();

        // Compare with mesh before save.
        CompareMeshes(&r_mesh, &r_mesh2);

        // Tidy up
        delete p_simulator;
        WntConcentration<2>::Destroy();
    }

    // A check that save and load works when a Voronoi tessellation is involved
    void TestMeshSurvivesSaveLoadWithVoronoiTessellation() throw (Exception)
    {
        // Create mesh
        double crypt_length = 22.0;
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, false);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.SetAreaBasedDampingConstant(true);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt, false, true);
        simulator.SetOutputDirectory("Crypt2DMeshArchive2");
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(LinearSpringWithVariableSpringConstantsForce<2>, p_variable_force);
        p_variable_force->SetEdgeBasedSpringConstant(true);
        simulator.AddForce(p_variable_force);

        // Run simulation
        simulator.Solve();

        // Save
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        // Load
        CryptSimulation2d* p_simulator;
        p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("Crypt2DMeshArchive2", 0.1);

        // Reset end time and run simulation
        p_simulator->SetEndTime(0.15);
        p_simulator->Solve();

        // Tidy up
        delete p_simulator;
        WntConcentration<2>::Destroy();
    }

    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        double crypt_length = cells_up*(sqrt(3)/2);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.SetOutputCellCyclePhases(true); // for coverage

        // We have a Wnt Gradient - but not Wnt dependent cells
        // so that the test runs quickly, but we test archiving of it!
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicStandardResult");
        simulator.SetEndTime(0.25);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25, which is also tested below
        // after a save and a load. (See #420, #479.)
        std::vector<double> node_28_location = simulator.GetNodeLocation(28);
        TS_ASSERT_DELTA(node_28_location[0], 3.7875, 1e-4);
        TS_ASSERT_DELTA(node_28_location[1], 0.0, 1e-4);

        std::vector<double> node_120_location = simulator.GetNodeLocation(120);
        TS_ASSERT_DELTA(node_120_location[0], 4.2035, 1e-4);
        TS_ASSERT_DELTA(node_120_location[1], 0.1033, 1e-4);

        // Test the Wnt concentration result
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(28)), 1.0, 1e-9);
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(120)), 0.9900, 1e-4);
        WntConcentration<2>::Destroy();

        // Check writing of Voronoi data
        OutputFileHandler handler("Crypt2DPeriodicStandardResult", false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellcyclephases.dat";

        NumericFileComparison comp(results_file, "crypt/test/data/CellCyclePhaseOutput/cellcyclephases.dat");
        TS_ASSERT(comp.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_file + " crypt/test/data/CellCyclePhaseOutput/cellcyclephases.dat").c_str()), 0);
    }

    // Testing Save
    void TestSave() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        double crypt_length = cells_up*(sqrt(3)/2);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        TS_ASSERT_EQUALS(simulator.GetIdentifier(), "CryptSimulation2d");

        simulator.SetOutputDirectory("Crypt2DPeriodicSaveAndLoad");

        // Our full end time is 0.25, here we run until 0.1 then load and run more below.
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception)
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        CryptSimulation2d* p_simulator1;

        // Make sure there is no existing WntConcentration instance before load
        WntConcentration<2>::Instance();
        WntConcentration<2>::Destroy();

        p_simulator1 = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("Crypt2DPeriodicSaveAndLoad", 0.1);
        p_simulator1->SetEndTime(0.2);
        p_simulator1->Solve();

        // Get mesh
        MutableMesh<2,2>& r_mesh1 = (static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator1->rGetCellPopulation())))->rGetMesh();

        // Save then reload, compare meshes either side
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator1);

        CryptSimulation2d* p_simulator2 = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("Crypt2DPeriodicSaveAndLoad", 0.2);
        MutableMesh<2,2>& r_mesh2 = (static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator2->rGetCellPopulation())))->rGetMesh();

        CompareMeshes(&r_mesh1, &r_mesh2);

        // Run a bit further...
        p_simulator2->SetEndTime(0.25);

        // Run simulation
        p_simulator2->Solve();

        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25 in the StandardResult test above.
        std::vector<double> node_28_location = p_simulator2->GetNodeLocation(28);
        TS_ASSERT_DELTA(node_28_location[0], 3.7875, 1e-4);
        TS_ASSERT_DELTA(node_28_location[1], 0.0, 1e-4);

        std::vector<double> node_120_location = p_simulator2->GetNodeLocation(120);
        TS_ASSERT_DELTA(node_120_location[0], 4.2035, 1e-4);
        TS_ASSERT_DELTA(node_120_location[1], 0.1033, 1e-4);

        // Test WntConcentration was set up correctly
        TS_ASSERT_EQUALS(WntConcentration<2>::Instance()->IsWntSetUp(), true);

        // Test the WntConcentration result
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(p_simulator2->rGetCellPopulation().GetCellUsingLocationIndex(28)), 1.0, 1e-9);
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(p_simulator2->rGetCellPopulation().GetCellUsingLocationIndex(120)), 0.9900, 1e-4);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
        WntConcentration<2>::Destroy();
    }

    /*
     * Cells are compressed and are trying to spread out. So the cells
     * at the bottom try to push across y=0 but are prevented from doing so.
     * This is not covered in previous tests because cells are all at age = 0
     * and try to come together to simulate birth.
     *
     * It is potentially an expensive test computationally, because all
     * Wnt cells have to run cell-cycle models for a large time
     * to be 'mature' cells which won't shrink together.
     * Limited this by using only four cells of minimum age.
     */
    void TestWntCellsCannotMoveAcrossYEqualsZero() throw (Exception)
    {
        // Create mesh
        double crypt_length = 22.0;
        unsigned cells_across = 2;
        unsigned cells_up = 2;
        double crypt_width = 0.5; // Make this bigger if want to visualise output
        unsigned thickness_of_ghost_layer = 1;

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Cover the write Voronoi data method
        crypt.SetOutputVoronoiData(true);
        crypt.SetOutputCellPopulationVolumes(true);
        crypt.SetOutputCellVolumes(true);
        crypt.SetOutputCellAncestors(true);
        crypt.SetOutputCellAges(true);

        AbstractCellPopulation<2>::Iterator cell_iterator = crypt.Begin();
        cell_iterator->SetBirthTime(-1.0);   // Make cell-cycle models do minimum work
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);

        boost::shared_ptr<AbstractCellProperty> p_apc1(crypt.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(crypt.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_label(crypt.GetCellPropertyRegistry()->Get<CellLabel>());

        cell_iterator->AddCellProperty(p_label);
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(p_apc1);
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(p_bcat1);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DWntMatureCells");
        simulator.SetEndTime(0.01);

        // Create a force law and pass it to the simulation
        MAKE_PTR(LinearSpringWithVariableSpringConstantsForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // If you want to visualize this use the 'notcylindrical' option
        // (it is too small for it to figure out what's happening on its own)

        // Cover exceptions
        TS_ASSERT_THROWS_THIS(simulator.rGetCellPopulation().GetCellMutationStateCount(),
                              "Call SetOutputCellMutationStates(true) before using this function");
        simulator.rGetCellPopulation().SetOutputCellMutationStates(true);

        TS_ASSERT_THROWS_THIS(simulator.rGetCellPopulation().rGetCellProliferativeTypeCount(),
                              "Call SetOutputCellProliferativeTypes(true) before using this function");
        simulator.rGetCellPopulation().SetOutputCellProliferativeTypes(true);

        TS_ASSERT_THROWS_THIS(simulator.rGetCellPopulation().rGetCellCyclePhaseCount(),
                              "Call SetOutputCellCyclePhases(true) before using this function");
        crypt.SetOutputCellCyclePhases(true);

        // Run simulation
        simulator.Solve();

        // Check that nothing has moved below y=0
        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            TS_ASSERT_LESS_THAN(-1e-15, crypt.GetLocationOfCellCentre(*cell_iter)[1]);
        }

        std::vector<unsigned> cell_mutation_state_count = simulator.rGetCellPopulation().GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_state_count.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[0], 2u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[2], 0u); // No APC two hit
        TS_ASSERT_EQUALS(cell_mutation_state_count[3], 1u);

        std::vector<unsigned> cell_type_count = simulator.rGetCellPopulation().rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count[0], 0u);
        TS_ASSERT_EQUALS(cell_type_count[1], 4u);
        TS_ASSERT_EQUALS(cell_type_count[2], 0u);

        // Check writing of Voronoi data
        OutputFileHandler handler("Crypt2DWntMatureCells", false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/voronoi.dat";

        NumericFileComparison comp(results_file, "crypt/test/data/Crypt2DWntMatureCells/voronoi.dat");
        TS_ASSERT(comp.CompareFiles(2e-6));

        // Tidy up
        WntConcentration<2>::Destroy();
    }


    void TestCellIdOutput() throw (Exception)
    {
        // Resetting the Maximum cell Id to zero (to account for previous tests)
        Cell::ResetMaxCellId();

        // Create mesh
        double crypt_length = 22.0;
        unsigned cells_across = 6;
        unsigned cells_up = 8;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);// true = mature cells

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-11.5);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Cover writing logged cell
        crypt.SetOutputCellIdData(true);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DCylindricalCellIdLogged");
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        // Check writing of cell data
        OutputFileHandler handler("Crypt2DCylindricalCellIdLogged", false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/loggedcell.dat";

        NumericFileComparison comp(results_file, "crypt/test/data/Crypt2DCylindricalCellIdLogged/loggedcell.dat");
        TS_ASSERT(comp.CompareFiles(2e-4));
    }

    // This is a strange test -- all cells divide within a quick time, it gives
    // good testing of the periodic boundaries though... [comment no longer valid?]
    void TestWithTysonNovakCells() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        double crypt_length = cells_up*(sqrt(3)/2);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<TysonNovakCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicTysonNovak");
        simulator.SetEndTime(0.05);
        simulator.SetDt(0.001);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Test that labelling a few cells doesn't make any difference to the simulation
        // and therefore log them in the visualizer files for the next test to check.
        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(ApcOneHitCellMutationState, p_apc1);
        MAKE_PTR(ApcTwoHitCellMutationState, p_apc2);
        MAKE_PTR(BetaCateninOneHitCellMutationState, p_bcat1);

        simulator.rGetCellPopulation().GetCellUsingLocationIndex(57)->AddCellProperty(p_label);
        simulator.rGetCellPopulation().GetCellUsingLocationIndex(56)->SetMutationState(p_apc1);
        simulator.rGetCellPopulation().GetCellUsingLocationIndex(51)->SetMutationState(p_apc2);
        simulator.rGetCellPopulation().GetCellUsingLocationIndex(63)->SetMutationState(p_bcat1);
        simulator.rGetCellPopulation().SetOutputCellMutationStates(true);

        // Run simulation
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_nodes,ghost_cells.size());
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 76u);
        TS_ASSERT_EQUALS(number_of_nodes, 124u);
    }

    void TestAddCellKiller() throw (Exception)
    {
        // Create mesh
        double crypt_length = 9.3;

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, &mesh, std::vector<unsigned>(), false, 0.0, 3.0, 6.5, 8.0);

        cells[60]->SetBirthTime(-50.0);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(mesh, cells);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        unsigned num_deaths = simulator.DoCellRemoval();
        unsigned num_births = simulator.DoCellBirth();

        TS_ASSERT_EQUALS(num_births, 1u);
        TS_ASSERT_EQUALS(num_deaths, 11u);
    }


    void TestCalculateCellDivisionVectorConfMesh() throw (Exception)
    {
        // Make a parent node
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);

        MutableMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);

        // Create cells
        std::vector<CellPtr> conf_cells;
        CryptCellsGenerator<TysonNovakCellCycleModel> cells_generator;
        cells_generator.Generate(conf_cells, &conf_mesh, std::vector<unsigned>(), true);

        // Create cell population
        MeshBasedCellPopulation<2> conf_crypt(conf_mesh, conf_cells);
        conf_crypt.SetMeinekeDivisionSeparation(0.1);

        AbstractCellPopulation<2>::Iterator conf_iter = conf_crypt.Begin();

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(conf_crypt);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeDivisionRestingSpringLength(0.9); // coverage
        simulator.AddForce(p_linear_force);

        c_vector<double, 2> daughter_location = simulator.CalculateCellDivisionVector(*conf_iter);
        c_vector<double, 2> new_parent_location = conf_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = conf_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), conf_crypt.GetMeinekeDivisionSeparation(), 1e-7);
    }


    void TestCalculateCellDivisionVectorConfMeshStemCell() throw (Exception)
    {
        // Make a parent node
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 0.0; // <- y=0
        Node<2>* p_node = new Node<2>(0u,location, false);
        MutableMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);

        // Create cells
        std::vector<CellPtr> conf_cells;
        CryptCellsGenerator<TysonNovakCellCycleModel> cells_generator;
        cells_generator.Generate(conf_cells, &conf_mesh, std::vector<unsigned>(), true);

        // Create cell population
        MeshBasedCellPopulation<2> conf_crypt(conf_mesh, conf_cells);

        AbstractCellPopulation<2>::Iterator conf_iter = conf_crypt.Begin();

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(conf_crypt);

        // Repeat two times for coverage
        // need vector from parent to daughter to have both +ve and -ve y component
        // different branches will execute to make sure daughter stays in crypt ie. +ve y component
        for (unsigned repetitions=0; repetitions<=1; repetitions++)
        {
            c_vector<double, 2> daughter_location = simulator.CalculateCellDivisionVector(*conf_iter);
            c_vector<double, 2> new_parent_location = conf_mesh.GetNode(0)->rGetLocation();
            c_vector<double, 2> parent_to_daughter = conf_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);

            // The parent stem cell should stay where it is and the daughter be introduced at positive y.

            TS_ASSERT_DELTA(new_parent_location[0], location[0], 1e-7);
            TS_ASSERT_DELTA(new_parent_location[1], location[1], 1e-7);
            TS_ASSERT_LESS_THAN_EQUALS(location[1], daughter_location[1]);
            TS_ASSERT_DELTA(norm_2(parent_to_daughter), conf_crypt.GetMeinekeDivisionSeparation(), 1e-7);
       }
    }


    void TestCalculateCellDivisionVectorCylindricalMesh() throw (Exception)
    {
        // Make a mesh
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);
        Cylindrical2dMesh cyl_mesh(6.0);
        cyl_mesh.AddNode(p_node);

        // Create cells
        std::vector<CellPtr> cyl_cells;
        CryptCellsGenerator<TysonNovakCellCycleModel> cells_generator;
        cells_generator.Generate(cyl_cells, &cyl_mesh, std::vector<unsigned>(), true);

        // Create cell population
        MeshBasedCellPopulation<2> cyl_crypt(cyl_mesh, cyl_cells);

        AbstractCellPopulation<2>::Iterator cyl_iter = cyl_crypt.Begin();

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(cyl_crypt);

        c_vector<double, 2> daughter_location = simulator.CalculateCellDivisionVector(*cyl_iter);
        c_vector<double, 2> new_parent_location = cyl_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = cyl_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), cyl_crypt.GetMeinekeDivisionSeparation(), 1e-7);
    }


    void TestCalculateCellDivisionVectorCylindricalMeshStemCell() throw (Exception)
    {
        // Make a mesh
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 0.0; // <- y=0
        Node<2>* p_node = new Node<2>(0u,location, false);
        Cylindrical2dMesh cyl_mesh(6.0);
        cyl_mesh.AddNode(p_node);

        // Create cells
        std::vector<CellPtr> cyl_cells;
        CryptCellsGenerator<TysonNovakCellCycleModel> cells_generator;
        cells_generator.Generate(cyl_cells, &cyl_mesh, std::vector<unsigned>(), true);

        // Create cell population
        MeshBasedCellPopulation<2> cyl_crypt(cyl_mesh, cyl_cells);

        AbstractCellPopulation<2>::Iterator cyl_iter = cyl_crypt.Begin();

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(cyl_crypt);

        c_vector<double,2> daughter_location = simulator.CalculateCellDivisionVector(*cyl_iter);
        c_vector<double,2> new_parent_location = cyl_mesh.GetNode(0)->rGetLocation();
        c_vector<double,2> parent_to_daughter = cyl_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);

        // The parent stem cell should stay where it is and the daughter be introduced at positive y.
        TS_ASSERT_DELTA(new_parent_location[0], location[0], 1e-7);
        TS_ASSERT_DELTA(new_parent_location[1], location[1], 1e-7);
        TS_ASSERT_LESS_THAN_EQUALS(location[1], daughter_location[1]);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), cyl_crypt.GetMeinekeDivisionSeparation(), 1e-7);
    }


    // Short test which sets mNoBirth for coverage
    void TestNoBirth() throw (Exception)
    {
        std::string output_directory = "Crypt2DCylindricalNoBirth";
        double crypt_length = 22.0;

        // Create mesh
        unsigned cells_across = 2;
        unsigned cells_up = 3;
        double crypt_width = 2.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);// true = mature cells

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(2.0); // long enough for a cell to be born were SetNoBirth not called

        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetUpdateCellPopulationRule(true);
        simulator.SetNoBirth(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_cells, cells_across*cells_up);
        TS_ASSERT_EQUALS(number_of_nodes, number_of_cells+thickness_of_ghost_layer*2*cells_across);
    }


    // Test death on a non-periodic mesh. Note that birth does occur too.
    void TestRandomDeathOnNonPeriodicCrypt() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 2;
        unsigned cells_up = 1;
        unsigned thickness_of_ghost_layer = 1;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DRandomDeathNonPeriodic");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_killer, (&crypt, 0.999996771));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // There should be no cells left at this time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 1u);
    }


    void TestUsingJiggledBottomSurface()
    {
        // Create mesh
        CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DJiggledBottomCells");
        simulator.SetEndTime(0.01);
        simulator.UseJiggledBottomCells();

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Move the first cell (which should be on y=0) down a bit
        AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
        TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0, 1e-6);

        // Move the cell (can't use the iterator for this as it is const)
        crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;
        TS_ASSERT_LESS_THAN(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0);

        // Run simulation
        simulator.Solve();

        // The cell should have been pulled up, but not above y=0. However it should
        // then been moved to above y=0 by the jiggling
        TS_ASSERT_LESS_THAN(0.0, crypt.GetLocationOfCellCentre(*cell_iter)[1]);
    }

    /**
     * Test that the cell count vectors are correctly initialized when a
     * simulation is saved then loaded.
     */
    void TestCellCountInitialization()
    {
        // Create mesh
        CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);
        TS_ASSERT_EQUALS(cells.size(), 16u);

        // Bestow mutations on some cells
        cells[0]->SetMutationState(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        cells[1]->SetMutationState(CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>());
        cells[2]->SetMutationState(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
        cells[3]->SetMutationState(CellPropertyRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>());
        cells[4]->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellLabel>());
        cells[2]->SetBirthTime(1.5-(cells[2]->GetCellCycleModel()->GetStemCellG1Duration()
                                   + cells[2]->GetCellCycleModel()->GetSG2MDuration()));

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        crypt.SetOutputCellMutationStates(true);
        crypt.SetOutputCellProliferativeTypes(true);
        crypt.SetOutputCellCyclePhases(true);

        // Each cell count has been initialized to the correct size and computed
        std::vector<unsigned> cell_mutation_state_count1 = crypt.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_state_count1.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_state_count1[0], 13u);
        TS_ASSERT_EQUALS(cell_mutation_state_count1[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count1[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count1[3], 1u);

        // However, using the cell mutation objects, they have counted!
        boost::shared_ptr<CellPropertyRegistry> p_registry = crypt.GetCellPropertyRegistry();
        TS_ASSERT_EQUALS(p_registry->Get<WildTypeCellMutationState>()->GetCellCount(), 13u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcTwoHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<BetaCateninOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<CellLabel>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApoptoticCellProperty>()->GetCellCount(), 0u);

        std::vector<unsigned> cell_type_count1 = crypt.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count1.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count1[0], 0u);
        TS_ASSERT_EQUALS(cell_type_count1[1], 0u);
        TS_ASSERT_EQUALS(cell_type_count1[2], 0u);

        std::vector<unsigned> cell_cycle_phase_count1 = crypt.rGetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phase_count1.size(), 5u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count1[0], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count1[1], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count1[2], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count1[3], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count1[4], 0u);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);

        simulator.SetOutputDirectory("TestMutationStateCellCount");
        simulator.SetEndTime(1.0);
        simulator.UseJiggledBottomCells();

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        // Each cell count has now been computed since WriteCellResultsToFiles() has been called
        std::vector<unsigned> cell_mutation_state_count3 = simulator.rGetCellPopulation().GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_state_count3.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_state_count3[0], 14u);
        TS_ASSERT_EQUALS(cell_mutation_state_count3[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count3[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count3[3], 1u);

        p_registry = simulator.rGetCellPopulation().GetCellPropertyRegistry();
        TS_ASSERT_EQUALS(p_registry->Get<WildTypeCellMutationState>()->GetCellCount(), 14u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcTwoHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<BetaCateninOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<CellLabel>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApoptoticCellProperty>()->GetCellCount(), 0u);

        std::vector<unsigned> cell_type_count3 = crypt.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count3.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count3[0], 4u);
        TS_ASSERT_EQUALS(cell_type_count3[1], 13u);
        TS_ASSERT_EQUALS(cell_type_count3[2], 0u);

        std::vector<unsigned> cell_cycle_phase_count3 = crypt.rGetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phase_count3.size(), 5u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count3[0], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count3[1], 1u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count3[2], 7u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count3[3], 7u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count3[4], 2u);

        // Save the simulation
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        // Check the original simulation's counts haven't changed
        p_registry = simulator.rGetCellPopulation().GetCellPropertyRegistry();
        TS_ASSERT_EQUALS(p_registry->Get<WildTypeCellMutationState>()->GetCellCount(), 14u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcTwoHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<BetaCateninOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApoptoticCellProperty>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry->Get<CellLabel>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApoptoticCellProperty>()->GetCellCount(), 0u);

        // Load the simulation
        CryptSimulation2d* p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("TestMutationStateCellCount", 1.0);

        // Check the original simulation's counts haven't changed
        p_registry = simulator.rGetCellPopulation().GetCellPropertyRegistry();
        TS_ASSERT_EQUALS(p_registry->Get<WildTypeCellMutationState>()->GetCellCount(), 14u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcTwoHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<BetaCateninOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<CellLabel>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApoptoticCellProperty>()->GetCellCount(), 0u);

        // In the loaded simulation, we want the various cell counts to be saved
        // (so that simulations which quit when a certain population is removed don't stop too soon)
        std::vector<unsigned> cell_mutation_state_count4 = p_simulator->rGetCellPopulation().GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_state_count4.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_state_count4[0], 14u);
        TS_ASSERT_EQUALS(cell_mutation_state_count4[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count4[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count4[3], 1u);

        p_registry = p_simulator->rGetCellPopulation().GetCellPropertyRegistry();
        TS_ASSERT_EQUALS(p_registry->Get<WildTypeCellMutationState>()->GetCellCount(), 14u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcTwoHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<BetaCateninOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<CellLabel>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApoptoticCellProperty>()->GetCellCount(), 0u);

        std::vector<unsigned> cell_type_count4 = crypt.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count4.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count4[0], 4u);
        TS_ASSERT_EQUALS(cell_type_count4[1], 13u);
        TS_ASSERT_EQUALS(cell_type_count4[2], 0u);

        std::vector<unsigned> cell_cycle_phase_count4 = crypt.rGetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phase_count4.size(), 5u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count4[0], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count4[1], 1u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count4[2], 7u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count4[3], 7u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count4[4], 2u);

        // Run simulation for a further time
        p_simulator->SetEndTime(2.0);
        p_simulator->Solve();

        // The cell mutation state count has now been computed since WriteCellResultsToFiles() has been called
        std::vector<unsigned> cell_mutation_state_count5 = p_simulator->rGetCellPopulation().GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_state_count5.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_state_count5[0], 16u);
        TS_ASSERT_EQUALS(cell_mutation_state_count5[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count5[2], 2u);
        TS_ASSERT_EQUALS(cell_mutation_state_count5[3], 1u);

        p_registry = p_simulator->rGetCellPopulation().GetCellPropertyRegistry();
        TS_ASSERT_EQUALS(p_registry->Get<WildTypeCellMutationState>()->GetCellCount(), 16u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApcTwoHitCellMutationState>()->GetCellCount(), 2u);
        TS_ASSERT_EQUALS(p_registry->Get<BetaCateninOneHitCellMutationState>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(crypt.GetCellPropertyRegistry()->Get<CellLabel>()->GetCellCount(), 1u);
        TS_ASSERT_EQUALS(p_registry->Get<ApoptoticCellProperty>()->GetCellCount(), 0u);

        std::vector<unsigned> cell_type_count5 = p_simulator->rGetCellPopulation().rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count5.size(), 3u);
        TS_ASSERT_EQUALS(cell_type_count5[0], 4u);
        TS_ASSERT_EQUALS(cell_type_count5[1], 16u);
        TS_ASSERT_EQUALS(cell_type_count5[2], 0u);

        std::vector<unsigned> cell_cycle_phase_count5 = p_simulator->rGetCellPopulation().rGetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phase_count5.size(), 5u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count5[0], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count5[1], 3u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count5[2], 5u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count5[3], 6u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count5[4], 6u);

        // Tidy up
        delete p_simulator;
    }

    void TestWriteBetaCatenin() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 5;
        unsigned cells_up = 4;
        unsigned thickness_of_ghost_layer = 1;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        double crypt_length = cells_up*(sqrt(3)/2);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisOne> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, false);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("CryptBetaCatenin");
        simulator.SetEndTime(0.01);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        // Check writing of beta-catenin data
        OutputFileHandler handler("CryptBetaCatenin", false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizbetacatenin";
        std::string results_setup_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizsetup";

        NumericFileComparison comp_bcat(results_file, "crypt/test/data/CryptBetaCatenin/results.vizbetacatenin");
        TS_ASSERT(comp_bcat.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_file + " crypt/test/data/CryptBetaCatenin/results.vizbetacatenin").c_str()), 0);

        TS_ASSERT_EQUALS(system(("diff " + results_setup_file + " crypt/test/data/CryptBetaCatenin/results.vizsetup").c_str()), 0);

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    void TestCryptSimulation2DParameterOutput() throw (Exception)
    {
        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);// true = mature cells

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.SetOutputCellMutationStates(true);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        ///\todo #1453 add an killer and test the output is correct
        std::string output_directory = "TestCryptSimulation2dOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);
        out_stream parameter_file = output_file_handler.OpenOutputFile("crypt_sim_2d_results.parameters");
        simulator.OutputSimulationParameters(parameter_file);
        parameter_file->close();

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "crypt_sim_2d_results.parameters  crypt/test/data/TestCryptSimulationOutputParameters/crypt_sim_2d_results.parameters").c_str()), 0);

        ///\todo check output of simulator.OutputSimulationSetup();
    }

    void TestAncestorCryptSimulations() throw (Exception)
    {
        std::string output_directory = "AncestorCrypt";
        double crypt_length = 22.0;

        // Create mesh
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> temp_cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(temp_cells, p_mesh, std::vector<unsigned>(), true, 0.3, 2.0, 3.0, 4.0, true);

        // This awkward way of setting up the cells is a result of #430
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            cells.push_back(temp_cells[location_indices[i]]);
        }

        // Set up crypt
        MeshBasedCellPopulationWithGhostNodes<2>* p_crypt = new MeshBasedCellPopulationWithGhostNodes<2>(*p_mesh, cells, location_indices, false, 30.0); // Last parameter adjusts Ghost spring stiffness in line with the linear_force later on

        // Set simulation to output cell types and cell ancestors
        p_crypt->SetOutputCellMutationStates(true);
        p_crypt->SetOutputCellAncestors(true);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(*p_crypt, false, false);
        simulator.SetOutputDirectory(output_directory);

        // Set length of simulation here
        double time_of_each_run = 10.0*simulator.GetDt(); // for each run
        simulator.SetEndTime(time_of_each_run);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(30.0); //normally 15.0;
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&(simulator.rGetCellPopulation()), crypt_length));
        simulator.AddCellKiller(p_killer);

        simulator.UseJiggledBottomCells();
        simulator.SetBottomCellAncestors();

        // Run simulation
        simulator.Solve();

        // ... and checking visualization of labelled cells against previous run
        OutputFileHandler handler("AncestorCrypt", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizancestors";
        std::string results_file3 = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";

        NumericFileComparison comp_nodes(results_file1, "crypt/test/data/AncestorCrypt/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles());
        NumericFileComparison comp_ans(results_file2, "crypt/test/data/AncestorCrypt/results.vizancestors");
        TS_ASSERT(comp_ans.CompareFiles());
        NumericFileComparison comp_celltypes(results_file3, "crypt/test/data/AncestorCrypt/results.vizcelltypes");
        TS_ASSERT(comp_celltypes.CompareFiles());

        TS_ASSERT_EQUALS(system(("diff " + results_file1 + " crypt/test/data/AncestorCrypt/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_file2 + " crypt/test/data/AncestorCrypt/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_file3 + " crypt/test/data/AncestorCrypt/results.vizcelltypes").c_str()), 0);

        // Tidy up
        WntConcentration<2>::Destroy();
        delete p_crypt;
    }

    void noTestGeneralisedLinearSpringForceWithSpringConstantsForIngeBCatCells()
    {
        // Create mesh
        double crypt_length = 1.1*12.0*sqrt(3)/2.0;
        CylindricalHoneycombMeshGenerator generator(6, 12, 0, 1.1);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, false);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells,location_indices);

        // Set up WntConcentration and associate with cell population
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create simulation
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("TestLinearSpringWithVariableSpringConstantsForceModified");
        simulator.SetEndTime(0.5);

        crypt.CreateVoronoiTessellation(); // this method is normally called in a simulation loop

        // Create force and associate with simulation
        MAKE_PTR(LinearSpringWithVariableSpringConstantsForce<2>, p_linear_force);
        p_linear_force->SetBetaCateninSprings(true);
        simulator.AddForce(p_linear_force);

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Tidy up
        WntConcentration<2>::Destroy();
        SimulationTime::Destroy();
    }
};

#endif /*TESTCRYPTSIMULATION2DWITHMESHBASEDCELLPOPULATION_HPP_*/