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

#ifndef TESTCRYPTSIMULATION2DBOUNDARYCONDITION_HPP_
#define TESTCRYPTSIMULATION2DBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CryptSimulation2dBoundaryCondition.hpp"
#include "WntConcentration.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class TestCryptSimulation2dBoundaryCondition : public AbstractCellBasedTestSuite
{
public:

	void TestSetAndGetUseJiggledBottomCells() throw (Exception)
	{
		CryptSimulation2dBoundaryCondition boundary_condition(NULL);
		TS_ASSERT_EQUALS(boundary_condition.GetUseJiggledBottomCells(), false);

		boundary_condition.SetUseJiggledBottomCells(true);
		TS_ASSERT_EQUALS(boundary_condition.GetUseJiggledBottomCells(), true);
	}

    void TestConstructorWithCellPopulation() throw (Exception)
    {
        // Create a mesh-based cell population
        CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Pass this cell population to a boundary condition object
        CryptSimulation2dBoundaryCondition boundary_condition(&crypt);

        // Test access to the cell population
        const AbstractCellPopulation<2>* p_population = boundary_condition.GetCellPopulation();
        TS_ASSERT(p_population != NULL);
    }

    void TestOutputParameters() throw(Exception)
    {
        std::string output_directory = "TestOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        CryptSimulation2dBoundaryCondition boundary_condition(NULL);
        TS_ASSERT_EQUALS(boundary_condition.GetIdentifier(), "CryptSimulation2dBoundaryCondition");

        out_stream boundary_condition_parameter_file = output_file_handler.OpenOutputFile("results.parameters");
        boundary_condition.OutputCellPopulationBoundaryConditionParameters(boundary_condition_parameter_file);
        boundary_condition_parameter_file->close();

        std::string boundary_condition_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + boundary_condition_results_dir + "results.parameters crypt/test/data/TestCryptSimulation2dBoundaryCondition/results.parameters").c_str()), 0);

        // Test OutputCellPopulationBoundaryConditionInfo() method
        out_stream boundary_condition_info_file = output_file_handler.OpenOutputFile("results.info");
        boundary_condition.OutputCellPopulationBoundaryConditionInfo(boundary_condition_info_file);
        boundary_condition_info_file->close();

        TS_ASSERT_EQUALS(system(("diff " + boundary_condition_results_dir + "results.info crypt/test/data/TestCryptSimulation2dBoundaryCondition/results.info").c_str()), 0);
    }

    void TestImposeBoundaryConditionWithNoWntOrJiggling() throw(Exception)
	{
		// Create a cell population
		CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		std::vector<CellPtr> cells;
		CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
		cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

		MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

		// Store the node locations
        std::vector<c_vector<double, 2> > node_locations_before(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
        	node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

		// Now move the first cell (which should be on y=0) down a bit
		AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
		TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0, 1e-6);
		crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;
		TS_ASSERT_LESS_THAN(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0);

		// Create a boundary condition object
		CryptSimulation2dBoundaryCondition boundary_condition(&crypt);
		boundary_condition.SetUseJiggledBottomCells(false);

		// Impose the boundary condition without a Wnt stimulus or jiggled bottom cells
        boundary_condition.ImposeBoundaryCondition(node_locations_before);

		// The first cell should have been moved back by the boundary condition object
		TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0, 1e-4);

		// The nodes should all now be at their original locations
		std::vector<c_vector<double, 2> > node_locations_after(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
        	node_locations_after[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
        	TS_ASSERT_DELTA(node_locations_before[node_index](0), node_locations_after[node_index](0), 1e-3);
        	TS_ASSERT_DELTA(node_locations_before[node_index](1), node_locations_after[node_index](1), 1e-3);
        }
	}

    void TestImposeBoundaryConditionWithNoWntButWithJiggling() throw(Exception)
	{
		// Create a cell population
		CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		std::vector<CellPtr> cells;
		CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
		cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

		MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

		// Store the node locations
        std::vector<c_vector<double, 2> > node_locations_before(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
        	node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        // Move the first cell (which should be on y=0, and is a stem cell) down a bit
		AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
		crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;

        // Also move the second cell (which should be on y=0, and we make a transit cell)
		++cell_iter;
		TS_ASSERT_EQUALS(cell_iter->GetCellCycleModel()->GetCellProliferativeType(), STEM);
		cell_iter->GetCellCycleModel()->SetCellProliferativeType(TRANSIT);
		TS_ASSERT_EQUALS(cell_iter->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
		TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0, 1e-6);
		crypt.GetNode(1)->rGetModifiableLocation()[1] = -0.1;
		TS_ASSERT_LESS_THAN(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0);

		// Create a boundary condition object
		CryptSimulation2dBoundaryCondition boundary_condition(&crypt);
		boundary_condition.SetUseJiggledBottomCells(true);

		// Impose the boundary condition without a Wnt stimulus but with jiggled bottom cells
        boundary_condition.ImposeBoundaryCondition(node_locations_before);

		/*
		 * The first cell should have been pulled up to y=0 since it is a stem cell and
		 * there is no Wnt stimulus. It should then be unaffected by the jiggling. The
		 * second cell should not have been pulled up since it is a stem cell, but should
		 * then have been moved to above y=0 by the jiggling.
		 */
        AbstractCellPopulation<2>::Iterator cell_iter2 = crypt.Begin();
		TS_ASSERT_DELTA(0.0, crypt.GetLocationOfCellCentre(*cell_iter2)[1], 1e-3);
		++cell_iter2;
		TS_ASSERT_LESS_THAN(0.0, crypt.GetLocationOfCellCentre(*cell_iter2)[1]);
	}

    void TestImposeBoundaryConditionWithWntButNoJiggling() throw(Exception)
	{
		// Create a cell population
		CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		std::vector<CellPtr> cells;
		CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
		cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

		MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

		// Store the node locations
        std::vector<c_vector<double, 2> > node_locations_before(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
        	node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

		// Move the first cell (which should be on y=0) down a bit
		AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
		crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;
		c_vector<double, 2> perturbed_location = crypt.GetLocationOfCellCentre(*cell_iter);

		// Create a boundary condition object
		CryptSimulation2dBoundaryCondition boundary_condition(&crypt);
		boundary_condition.SetUseJiggledBottomCells(false);

		// Set up a WntConcentration singleton object
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt.GetWidth(1));
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

		// Impose the boundary condition with a Wnt stimulus but without jiggled bottom cells
        boundary_condition.ImposeBoundaryCondition(node_locations_before);

		/*
		 * The first cell should not have been moved by the boundary condition object,
		 * and hence should remain in its perturbed location.
		 */
		c_vector<double, 2> location_after = crypt.GetLocationOfCellCentre(*cell_iter);
		TS_ASSERT_DELTA(location_after[0], perturbed_location[0], 1e-3);
		TS_ASSERT_DELTA(location_after[1], location_after[1], 1e-3);


		// Tidy up
		WntConcentration<2>::Destroy();
	}

    void TestVerifyBoundaryCondition() throw(Exception)
	{
		// Create a cell population
		CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		std::vector<CellPtr> cells;
		CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
		cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

		MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

		// Store the node locations
        std::vector<c_vector<double, 2> > node_locations_before(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
        	node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

		// Now move the first cell (which should be on y=0) down a bit
		AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
		crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;

		// Create a boundary condition object
		CryptSimulation2dBoundaryCondition boundary_condition(&crypt);

		// Before imposing the boundary condition, it should not be satisfied
		TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), false);

		boundary_condition.ImposeBoundaryCondition(node_locations_before);

		// After imposing the boundary condition, it should be satisfied
		TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), true);
	}

    void TestArchiving() throw (Exception)
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false); // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CryptSimulation2dBoundaryCondition.arch";

        {
            // Create an output archive
        	CryptSimulation2dBoundaryCondition boundary_condition(NULL);
    		boundary_condition.SetUseJiggledBottomCells(true);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            AbstractCellPopulationBoundaryCondition<2>* const p_boundary_condition = &boundary_condition;
            output_arch << p_boundary_condition;
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCellPopulationBoundaryCondition<2>* p_boundary_condition;

            // Restore from the archive
            input_arch >> p_boundary_condition;

            // Test we have restored the plane geometry correctly
            TS_ASSERT_EQUALS(static_cast<CryptSimulation2dBoundaryCondition*>(p_boundary_condition)->GetUseJiggledBottomCells(), true);

            // Tidy up
            delete p_boundary_condition;
       }
    }
};

#endif /* TESTCRYPTSIMULATION2DBOUNDARYCONDITION_HPP_ */
