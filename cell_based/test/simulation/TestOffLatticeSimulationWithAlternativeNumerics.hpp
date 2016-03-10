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

#ifndef TESTOFFLATTICESIMULATIONWITHALTERNATIVENUMERICS_HPP_
#define TESTOFFLATTICESIMULATIONWITHALTERNATIVENUMERICS_HPP_

#include <cxxtest/TestSuite.h>

#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "RK4NumericalMethod.hpp"

/**
* Test simulations for the new numerical methods
*/

class TestOffLatticeSimulationWithAlternativeNumerics : public AbstractCellBasedTestSuite
{

public:

    void TestDefaultNodeBasedSimulationWithArchiving() throw (Exception)
    {
    	std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false,  0.25, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, -0.25, 0.0, 0.0));
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
        
        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);
	
       	// Make a simulation with the default numerical method
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("NodeBasedFEWithArchiving");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(1.0);

        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Run and archive the simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        CellBasedSimulationArchiver<3,OffLatticeSimulation<3> >::Save(&simulator);

        // Load and continue the simulation
        OffLatticeSimulation<3>* p_new_simulator;
        p_new_simulator = CellBasedSimulationArchiver<3,OffLatticeSimulation<3> >::Load("NodeBasedFEWithArchiving", 1.0);
        p_new_simulator->SetEndTime(2.0);
        TS_ASSERT_THROWS_NOTHING(p_new_simulator->Solve());
    }


    void TestRK4NodeBasedSimulationWithArchiving() throw (Exception)
    {
    	std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false,  0.25, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, -0.25, 0.0, 0.0));
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
        
        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);
	
       	// Make a simulation with the RK4 numerical method
       	MAKE_PTR(RK4NumericalMethod<3>, rk4Method);
        OffLatticeSimulation<3> simulator(node_based_cell_population, false, true, rk4Method);
        simulator.SetOutputDirectory("NodeBasedRK4WithArchiving");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(1.0);

        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Run and archive the simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        CellBasedSimulationArchiver<3,OffLatticeSimulation<3> >::Save(&simulator);

        // Load and continue the simulation
        OffLatticeSimulation<3>* p_new_simulator;
        p_new_simulator = CellBasedSimulationArchiver<3,OffLatticeSimulation<3> >::Load("NodeBasedRK4WithArchiving", 1.0);
        p_new_simulator->SetEndTime(2.0);
        TS_ASSERT_THROWS_NOTHING(p_new_simulator->Solve());
    }
};
  
#endif /*TESTOFFLATTICESIMULATIONWITHALTERNATIVENUMERICS_HPP_*/