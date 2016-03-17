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
#include "CellBasedEventHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "OffLatticeSimulation.hpp"
#include "Warnings.hpp"

#include "AbstractOffLatticeCellPopulation.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DiffusionForce.hpp"
#include "StochasticDurationCellCycleModel.hpp"

#include "ForwardEulerNumericalMethod.hpp"
#include "RK4NumericalMethod.hpp"

#include "PetscSetupAndFinalize.hpp"


/**
* Test simulations for the new numerical methods. 
* Each test does a short simulation, then archives, then forces a stepsize error 
*/

class TestOffLatticeSimulationWithAlternativeNumerics : public AbstractCellBasedTestSuite
{

public:

    void TestDefaultNodeBasedSimulationWithParticles() //throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Make a typical node based simulation with particles
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0,  true,  0.25, 0.0, 0.0));
        nodes.push_back(new Node<3>(1,  true, -0.25, 0.0, 0.0));
        nodes.push_back(new Node<3>(2,  false, 0.0, 0.25, 0.0));
        nodes.push_back(new Node<3>(3,  false, 0.0, -0.25, 0.0));

        MAKE_PTR(NodesOnlyMesh<3>, p_mesh);
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        // Specify which nodes correspond to cells (the rest are particles)
        std::vector<unsigned> location_indices;
        location_indices.push_back(0u);
        location_indices.push_back(1u);
        location_indices.push_back(2u);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);       
        NodeBasedCellPopulationWithParticles<3> cell_population(*p_mesh, cells, location_indices);

       	// Make a simulation using the default numerical method (FE)
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedFE");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(0.1);

        // Add a spring force
        MAKE_PTR(DiffusionForce<3>, p_diff_force);
        simulator.AddForce(p_diff_force);

        // Run and archive the simulation
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        CellBasedSimulationArchiver<3,OffLatticeSimulation<3> >::Save(&simulator);

        // Load the simulation and exceed the Absolute Movement Threshold
        OffLatticeSimulation<3>* p_new_simulator;
        p_new_simulator = CellBasedSimulationArchiver<3,OffLatticeSimulation<3> >::Load("NodeBasedFE", 0.1);
        AbstractCellPopulation<3,3>* p_population = &(p_new_simulator->rGetCellPopulation());
        dynamic_cast<AbstractOffLatticeCellPopulation<3,3>*>(p_population)->SetAbsoluteMovementThreshold(0.001);
        p_new_simulator->SetEndTime(0.2);
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_CONTAINS(p_new_simulator->Solve(), "which is more than the AbsoluteMovementThreshold:");

        for (unsigned i=0; i<nodes.size();i++)
        {
            delete nodes[i];
        }
    }


    void TestRK4NodeBasedSimulation() //throw (Exception)
    {   
        EXIT_IF_PARALLEL;

        // Make a typical node based simulation with particles
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0,  true,  0.25, 0.0, 0.0));
        nodes.push_back(new Node<3>(1,  true, -0.25, 0.0, 0.0));
        nodes.push_back(new Node<3>(2,  false, 0.0, 0.25, 0.0));
        nodes.push_back(new Node<3>(3,  false, 0.0, -0.25, 0.0));

        MAKE_PTR(NodesOnlyMesh<3>, p_mesh);
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        // Specify which nodes correspond to cells (the rest are particles)
        std::vector<unsigned> location_indices;
        location_indices.push_back(0u);
        location_indices.push_back(1u);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);       
        NodeBasedCellPopulationWithParticles<3> cell_population(*p_mesh, cells, location_indices);

        // Make a simulation using the RK4 numerical method
        MAKE_PTR(RK4NumericalMethod<3>, rk4Method);
        OffLatticeSimulation<3> simulator(cell_population, false, true, rk4Method);
        simulator.SetOutputDirectory("NodeBasedRK");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(0.1);

        // Add a spring force
        MAKE_PTR(DiffusionForce<3>, p_diff_force);
        simulator.AddForce(p_diff_force);

        // Run and archive the simulation
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        CellBasedSimulationArchiver<3,OffLatticeSimulation<3> >::Save(&simulator);

        // Load the simulation and exceed the Absolute Movement Threshold
        OffLatticeSimulation<3>* p_new_simulator;
        p_new_simulator = CellBasedSimulationArchiver<3,OffLatticeSimulation<3> >::Load("NodeBasedRK", 0.1);
        AbstractCellPopulation<3,3>* p_population = &(p_new_simulator->rGetCellPopulation());
        dynamic_cast<AbstractOffLatticeCellPopulation<3,3>*>(p_population)->SetAbsoluteMovementThreshold(0.001);
        p_new_simulator->SetEndTime(0.2);
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_CONTAINS(p_new_simulator->Solve(), "which is more than the AbsoluteMovementThreshold:"); 

        for (unsigned i=0; i<nodes.size();i++)
        {
            delete nodes[i];
        }
    }

    void TestDefaultVertexBasedSimulation() throw (Exception){
        
    	// Make a typical vertex based simulation
        HoneycombVertexMeshGenerator generator(5, 5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up a simulation with the default FE numerical method
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VertexBasedFE");
        simulator.SetSamplingTimestepMultiple(500);
        simulator.SetDt(1.0/500.0);
        simulator.SetEndTime(0.1);

        // Add a Nagai-Honda force
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run a simulation. Check that the appropriate step size warning is there (this one always seems to trigger???) 
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
        
        // Test archiving
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2> >::Save(&simulator);
        OffLatticeSimulation<2>* p_new_simulator;
        p_new_simulator = CellBasedSimulationArchiver<2,OffLatticeSimulation<2> >::Load("VertexBasedFE", 0.1);
        p_new_simulator->SetEndTime(0.2);
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_NOTHING(p_new_simulator->Solve());
        Warnings::QuietDestroy();
    }

    void TestRK4VertexBasedSimulation() throw (Exception){

        // Make a typical vertex based simulation
        HoneycombVertexMeshGenerator generator(5, 5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up a simulation with an RK4 numerical method
        MAKE_PTR(RK4NumericalMethod<2>, rk4Method);
        OffLatticeSimulation<2> simulator(cell_population, false, true, rk4Method);
        simulator.SetOutputDirectory("VertexBasedRK4");
        simulator.SetSamplingTimestepMultiple(500);
        simulator.SetDt(1.0/500.0);
        simulator.SetEndTime(0.1);

        // Add a Nagai-Honda force
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run a simulation. Check that the appropriate step size warning is there (this one always seems to trigger???) 
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
        
        // Test archiving
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2> >::Save(&simulator);
        OffLatticeSimulation<2>* p_new_simulator;
        p_new_simulator = CellBasedSimulationArchiver<2,OffLatticeSimulation<2> >::Load("VertexBasedRK4", 0.1);
        p_new_simulator->SetEndTime(0.2);
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_NOTHING(p_new_simulator->Solve());
        Warnings::QuietDestroy();
    }

    void TestDefaultMeshBasedWithGhostNodesSimulation() throw (Exception){
        
    	// Create a typical MeshBasedSimulationWithGhostNodes
		HoneycombMeshGenerator generator(6, 6, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        MeshBasedCellPopulationWithGhostNodes<2> cellPopulation(*p_mesh, cells, location_indices);
        cellPopulation.CreateVoronoiTessellation();
        
        // Create a simulation with the default FE numerical method
        OffLatticeSimulation<2> simulator(cellPopulation);
        simulator.SetOutputDirectory("MeshBasedWithGhostNodesFE");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetDt(1.0/120.0);
        simulator.SetEndTime(0.1);

        // Add a spring force law
        MAKE_PTR(GeneralisedLinearSpringForce<2>, pLinearForce);
        pLinearForce->SetCutOffLength(1.5);
        simulator.AddForce(pLinearForce);

        // Solve and check for errors
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test archiving and trigger an AbsoluteMovementThreshold error
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2> >::Save(&simulator);
        OffLatticeSimulation<2>* p_new_simulator;
        p_new_simulator = CellBasedSimulationArchiver<2,OffLatticeSimulation<2> >::Load("MeshBasedWithGhostNodesFE", 0.1);
        AbstractCellPopulation<2,2>* p_population = &(p_new_simulator->rGetCellPopulation());
        dynamic_cast<AbstractOffLatticeCellPopulation<2,2>*>(p_population)->SetAbsoluteMovementThreshold(0.001);
        p_new_simulator->SetEndTime(0.2);
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_CONTAINS(p_new_simulator->Solve(), "which is more than the AbsoluteMovementThreshold:");
    }

    void TestRK4MeshBasedWithGhostNodesSimulation() throw (Exception){

        // Create a typical MeshBasedSimulationWithGhostNodes
        HoneycombMeshGenerator generator(6, 6, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        MeshBasedCellPopulationWithGhostNodes<2> cellPopulation(*p_mesh, cells, location_indices);
        cellPopulation.CreateVoronoiTessellation();
        
        // Create a simulation with the default FE numerical method
        MAKE_PTR(RK4NumericalMethod<2>, rk4Method);
        OffLatticeSimulation<2> simulator(cellPopulation, false, true, rk4Method);
        simulator.SetOutputDirectory("MeshBasedWithGhostNodesRK4");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetDt(1.0/120.0);
        simulator.SetEndTime(0.1);

        // Add a spring force law
        MAKE_PTR(GeneralisedLinearSpringForce<2>, pLinearForce);
        pLinearForce->SetCutOffLength(1.5);
        simulator.AddForce(pLinearForce);

        // Solve and check for errors
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test archiving and trigger an AbsoluteMovementThreshold error
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2> >::Save(&simulator);
        OffLatticeSimulation<2>* p_new_simulator;
        p_new_simulator = CellBasedSimulationArchiver<2,OffLatticeSimulation<2> >::Load("MeshBasedWithGhostNodesRK4", 0.1);
        AbstractCellPopulation<2,2>* p_population = &(p_new_simulator->rGetCellPopulation());
        dynamic_cast<AbstractOffLatticeCellPopulation<2,2>*>(p_population)->SetAbsoluteMovementThreshold(0.001);
        p_new_simulator->SetEndTime(0.2);
        CellBasedEventHandler::Reset();
        TS_ASSERT_THROWS_CONTAINS(p_new_simulator->Solve(), "which is more than the AbsoluteMovementThreshold:");
    }
};
  
#endif /*TESTOFFLATTICESIMULATIONWITHALTERNATIVENUMERICS_HPP_*/