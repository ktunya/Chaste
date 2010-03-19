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
#ifndef TESTCELLSGENERATOR_HPP_
#define TESTCELLSGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "HoneycombMeshGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"

class TestCellsGenerator : public AbstractCellBasedTestSuite
{
public:

    void TestGenerateBasicWithFixedDurationGenerationBasedCellCycleModel() throw(Exception)
	{
		// Create mesh
		TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
		TetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);

		// Create cells
		std::vector<TissueCell> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

		// Test that cells were generated correctly
		TS_ASSERT_EQUALS(cells.size(), mesh.GetNumNodes());

		for (unsigned i=0; i<cells.size(); i++)
		{
			TS_ASSERT_DELTA(cells[i].GetBirthTime(), -(double)(i), 1e-9);
			TS_ASSERT_EQUALS(cells[i].GetCellCycleModel()->GetDimension(), 2u);
		}

		// Test with extra input argument
		std::vector<unsigned> location_indices;
		location_indices.push_back(2);
        location_indices.push_back(7);
        location_indices.push_back(9);

		std::vector<TissueCell> cells2;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, 3, location_indices);

        TS_ASSERT_EQUALS(cells2.size(), 3u);
        TS_ASSERT_DELTA(cells2[0].GetBirthTime(), -2.0, 1e-4);
        TS_ASSERT_DELTA(cells2[1].GetBirthTime(), -7.0, 1e-4);
        TS_ASSERT_DELTA(cells2[2].GetBirthTime(), -9.0, 1e-4);
	}

    void TestGenerateGivenLocationIndicesWithFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Use a mesh generator to generate some location indices corresponding to real cells
        HoneycombMeshGenerator mesh_generator(6, 7, 2, false);
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells again with basic
  		std::vector<TissueCell> cells;
  		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
  		TS_ASSERT_THROWS_THIS(cells_generator.GenerateBasic(cells, 83511u, location_indices),
  				"The size of the locationIndices vector must match the required number of output cells");
  		cells_generator.GenerateGivenLocationIndices(cells, location_indices);

  		// Test that cells were generated correctly
  		for (unsigned i=0; i<cells.size(); i++)
  		{
  			TS_ASSERT_DELTA(cells[i].GetBirthTime(), -(double)(location_indices[i]), 1e-9);
			TS_ASSERT_EQUALS(cells[i].GetCellCycleModel()->GetDimension(), 2u);
  		}
    }

    void TestCryptCellsGeneratorWithFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<TissueCell> cells;

        double y0 = 0.2;
        double y1 = 1.0;
        double y2 = 2.0;
        double y3 = 3.0;

        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, true, y0, y1, y2, y3);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        // Test that cells were generated correctly
        for (unsigned i=0; i<cells.size(); i++)
        {
            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = static_cast<FixedDurationGenerationBasedCellCycleModel*>(cells[i].GetCellCycleModel())->GetGeneration();

			TS_ASSERT_EQUALS(cells[i].GetCellCycleModel()->GetDimension(), 2u);

			if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
    			TS_ASSERT_EQUALS(cells[i].GetCellProliferativeType(), STEM);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
    			TS_ASSERT_EQUALS(cells[i].GetCellProliferativeType(), TRANSIT);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
    			TS_ASSERT_EQUALS(cells[i].GetCellProliferativeType(), TRANSIT);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
    			TS_ASSERT_EQUALS(cells[i].GetCellProliferativeType(), TRANSIT);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
    			TS_ASSERT_EQUALS(cells[i].GetCellProliferativeType(), DIFFERENTIATED);
            }
        }
    }

    void TestCryptCellsGeneratorWithStochasticDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<TissueCell> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, false);

        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        double y0 = 0.3;
        double y1 = 2.0;
        double y2 = 3.0;
        double y3 = 4.0;

        // Test that cells were generated correctly
        for (unsigned i=0; i<cells.size(); i++)
        {
            double height = p_mesh->GetNode(i)->rGetLocation()[1];
            unsigned generation = static_cast<StochasticDurationGenerationBasedCellCycleModel*>(cells[i].GetCellCycleModel())->GetGeneration();

            if (height <= y0)
            {
                TS_ASSERT_EQUALS(generation, 0u);
            }
            else if (height < y1)
            {
                TS_ASSERT_EQUALS(generation, 1u);
            }
            else if (height < y2)
            {
                TS_ASSERT_EQUALS(generation, 2u);
            }
            else if (height < y3)
            {
                TS_ASSERT_EQUALS(generation, 3u);
            }
            else
            {
                TS_ASSERT_EQUALS(generation, 4u);
            }

            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }

        // Create cells again with basic
		std::vector<TissueCell> new_cells;
		generator.GenerateBasic(new_cells, p_mesh->GetNumNodes());
		// Test that cells were generated correctly
		TS_ASSERT_EQUALS(new_cells.size(), p_mesh->GetNumNodes());

		for (unsigned i=0; i<new_cells.size(); i++)
		{
			TS_ASSERT_DELTA(new_cells[i].GetBirthTime(), -(double)(i), 1e-9);
		}
    }

    void TestTestCryptCellsGeneratorWithTysonNovakCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<TissueCell> cells;
        CryptCellsGenerator<TysonNovakCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, true);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        // Create cells again with basic
		std::vector<TissueCell> new_cells;
		generator.GenerateBasic(new_cells, p_mesh->GetNumNodes());

		// Test that cells were generated correctly
		TS_ASSERT_EQUALS(new_cells.size(), p_mesh->GetNumNodes());

		for (unsigned i=0; i<new_cells.size(); i++)
		{
			TS_ASSERT_DELTA(new_cells[i].GetBirthTime(), -(double)(i), 1e-9);
		}
    }


    void TestTestCryptCellsGeneratorWithWntCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<TissueCell> cells;
        CryptCellsGenerator<WntCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, false);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }

        // Create cells again with basic
		std::vector<TissueCell> new_cells;
		generator.GenerateBasic(new_cells, p_mesh->GetNumNodes());
		// Test that cells were generated correctly
		TS_ASSERT_EQUALS(new_cells.size(), p_mesh->GetNumNodes());

		for (unsigned i=0; i<new_cells.size(); i++)
		{
			TS_ASSERT_DELTA(new_cells[i].GetBirthTime(), -(double)(i), 1e-9);
		}
    }


    void TestTestCryptCellsGeneratorWithSimpleWntCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<TissueCell> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, false);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }
    }

    void TestTestCryptCellsGeneratorWithStochasticWntCellCycleModel() throw(Exception)
    {
        // Create mesh
        HoneycombMeshGenerator mesh_generator(5, 10, 0, false);
        TetrahedralMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        // Create cells
        std::vector<TissueCell> cells;
        CryptCellsGenerator<StochasticWntCellCycleModel> generator;
        generator.Generate(cells, p_mesh, location_indices, false);

        // Test that cells were generated correctly
        TS_ASSERT_EQUALS(cells.size(), p_mesh->GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_DELTA(cells[i].GetBirthTime(), 0.0, 1e-9);
        }

        // Create cells again with basic
		std::vector<TissueCell> new_cells;
		generator.GenerateBasic(new_cells, p_mesh->GetNumNodes());
		// Test that cells were generated correctly
		TS_ASSERT_EQUALS(new_cells.size(), p_mesh->GetNumNodes());

		for (unsigned i=0; i<new_cells.size(); i++)
		{
			TS_ASSERT_DELTA(new_cells[i].GetBirthTime(), -(double)(i), 1e-9);
		}
    }
    


    void TestCanDifferentiate() throw(Exception)
    {
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> gen1;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> gen2;
        CryptCellsGenerator<TysonNovakCellCycleModel> gen3;
        CryptCellsGenerator<WntCellCycleModel> gen4;
        CryptCellsGenerator<SimpleWntCellCycleModel> gen5;
        CryptCellsGenerator<StochasticWntCellCycleModel> gen6;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisOne> gen7;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo> gen8;
        
        TS_ASSERT_EQUALS(gen1.CanCellsTerminallyDifferentiate(), true);
        TS_ASSERT_EQUALS(gen2.CanCellsTerminallyDifferentiate(), true);
        TS_ASSERT_EQUALS(gen3.CanCellsTerminallyDifferentiate(), false);
        TS_ASSERT_EQUALS(gen4.CanCellsTerminallyDifferentiate(), false);
        TS_ASSERT_EQUALS(gen5.CanCellsTerminallyDifferentiate(), false);
        TS_ASSERT_EQUALS(gen6.CanCellsTerminallyDifferentiate(), false);
        TS_ASSERT_EQUALS(gen7.CanCellsTerminallyDifferentiate(), false);
        TS_ASSERT_EQUALS(gen8.CanCellsTerminallyDifferentiate(), false);
        
        CryptCellsGenerator<Alarcon2004OxygenBasedCellCycleModel> gen9;
        TS_ASSERT_THROWS_THIS( gen9.CanCellsTerminallyDifferentiate(), "Using an invalid cell cycle model for crypt simulations");
    }
};

#endif /*TESTCELLSGENERATOR_HPP_*/
