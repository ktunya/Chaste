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

#ifndef TESTLOGICCELLPOPULATION_HPP
#define TESTLOGICCELLPOPULATION_HPP

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "OffLatticeSimulation.hpp"
#include "RepulsionForce.hpp"

#include "LogicCell.hpp"
#include "EnvironmentalSignalTypes.hpp"
#include "DummyMutationState.hpp"
#include "CellCycleTestWithDivisionMechanism.hpp"
#include "DifferentiateIfSignalAboveThresh.hpp"
#include "PositionSignal.hpp"
#include "CellGrowth.hpp"


class TestLogicCellPopulation : public AbstractCellBasedTestSuite {

private:

   NodesOnlyMesh<3>* make3DCubeMesh(int cubeSide){

      std::vector< Node<3>* > nodes;
      int nodeIndex = 0;
      for (int i = 0; i < cubeSide; i++){
         for (int j = 0; j < cubeSide; j++){
            for (int k = 0; k < cubeSide; k++){
               Node<3>* newNode;
               newNode = new Node<3>(nodeIndex, false, k, j, i);
               nodes.push_back(newNode);
               nodeIndex++;
            }
         }
      }
      NodesOnlyMesh<3>* mesh = new NodesOnlyMesh<3>;
      mesh->ConstructNodesWithoutMesh(nodes, 5);

      for (unsigned i = 0; i<nodes.size(); i++)
      {
         delete nodes[i];
      }

      return mesh;
   }


   std::vector<boost::shared_ptr<PlaneBoundaryCondition<3> > > GetBoxWalls(NodeBasedCellPopulation<3>* cellPopulation, int sideOfInitialCellCube){

      std::vector<boost::shared_ptr<PlaneBoundaryCondition<3> > > walls;

      c_vector<double,3> point;
      c_vector<double,3> normal;

      point[0] = 0; point[1] = 0; point[2] = 0;
      normal[0] = 0; normal[1] = 0; normal[2] = -1;
      MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_base, (cellPopulation, point, normal));
      walls.push_back(p_base);
      point[0] = sideOfInitialCellCube*4; point[1] = 0; point[2] = 0;
      normal[0] = 1; normal[1] = 0; normal[2] = 0;
      MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_wall1, (cellPopulation, point, normal));
      walls.push_back(p_wall1);
      point[0] = sideOfInitialCellCube*4; point[1] = 0; point[2] = 0;
      normal[0] = 0; normal[1] = -1; normal[2] = 0;
      MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_wall2, (cellPopulation, point, normal));
      walls.push_back(p_wall2);
      point[0] = 0; point[1] = sideOfInitialCellCube*4; point[2] = 0;
      normal[0] = 0; normal[1] = 1; normal[2] = 0;
      MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_wall3, (cellPopulation, point, normal));
      walls.push_back(p_wall3);
      point[0] = 0; point[1] = sideOfInitialCellCube*4; point[2] = 0;
      normal[0] = -1; normal[1] = 0; normal[2] = 0;
      MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_wall4, (cellPopulation, point, normal));
      walls.push_back(p_wall4);

      return walls;
   }


   boost::shared_ptr<PlaneBasedCellKiller<3> > GetTopKiller(NodeBasedCellPopulation<3>* cellPopulation, int sideOfInitialCellCube){

      c_vector<double,3> point;
      c_vector<double,3> normal;

      point[0] = 0; point[1] = 0; point[2] = 10*sideOfInitialCellCube;
      normal[0] = 0; normal[1] = 0; normal[2] = 1;
      MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_topSloughing, (cellPopulation, point, normal));
      return p_topSloughing;
   }


public:

   void TestCellCount(){

      // Set up a cube of 8 cells
      int sideOfInitialCellCube = 2;
      int nCells = sideOfInitialCellCube * sideOfInitialCellCube * sideOfInitialCellCube;
      NodesOnlyMesh<3>* mesh  = make3DCubeMesh(sideOfInitialCellCube);

      std::vector<CellPtr> cells;
      cells.clear();
      cells.reserve(nCells);

      // Give each cell a test cell cycle logic, with fixed phase durations.
      for (int i=0; i < nCells; i++)
      {
         MAKE_PTR(DummyMutationState, pDummyMut);
         MAKE_PTR_ARGS(LogicCell, new_cell, (pDummyMut, false));

         CellCycleTestWithDivisionMechanism* cycleLogic = new CellCycleTestWithDivisionMechanism(new_cell.get(), CellCycle::G1, 0.0, 1.0);
         new_cell->SetLogic<CellCycle>(cycleLogic);

         cells.push_back(new_cell);
      }

      NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(*mesh, cells);
      OffLatticeSimulation<3>* simulation = new OffLatticeSimulation<3>(*cellPopulation);

      // Add a cell repulsion force 
      MAKE_PTR(RepulsionForce<3>, p_force);
      simulation->AddForce(p_force);

      // Set up and run simulation
      simulation->SetOutputDirectory("TestLogicCellPopulationCount");
      simulation->SetSamplingTimestepMultiple(120);
      simulation->SetDt(1.0/120.0);
      simulation->SetEndTime(12);
      simulation->Solve();     

      // Check that cells have divided exactly 3 times (64=8*(2^3))
      TS_ASSERT_EQUALS(simulation->rGetCellPopulation().GetNumNodes(), 64);
   }

   void TestColumnOfCells(){

      // Set up an initial cube of cells with side 3
      int sideOfInitialCellCube = 3;
      int nCells = sideOfInitialCellCube * sideOfInitialCellCube * sideOfInitialCellCube;
      NodesOnlyMesh<3>* mesh  = make3DCubeMesh(sideOfInitialCellCube);

      // Create cells
      std::vector<CellPtr> cells;
      cells.clear();
      cells.reserve(nCells);

      for (int i=0; i < nCells; i++)
      {

         MAKE_PTR(DummyMutationState, pDummyMut);
         MAKE_PTR_ARGS(LogicCell, new_cell, (pDummyMut, false));

         // Simple cell cycle model. All phase lengths = 2 hours.
         CellCycleTestWithDivisionMechanism* cycleLogic = new CellCycleTestWithDivisionMechanism(new_cell.get(), CellCycle::G1, 0.0, 2.0);
         new_cell->SetLogic<CellCycle>(cycleLogic);

         // A demo of cell growth logic. Cells halve in volume on division, grow at 0.5 microns per hour
         // up to a max radius of 2 microns, and begin with radius 1 micron  
         CellGrowth* growthLogic = new CellGrowth(new_cell.get(), Growth::GROWING, 1.0, 0.5, 2.0);
         new_cell->SetLogic<Growth>(growthLogic);

         // Add a simple distance based differentiation logic. S->S+S, TA->TA+TA and differentiation happens
         // 20 microns from the niche
         int startingState;
         if(i < sideOfInitialCellCube * sideOfInitialCellCube){
            // Bottom layer of cells is initialized to stem
            startingState = Differentiation::S;
         }else{
            startingState = Differentiation::TA;
         }
         DifferentiateIfSignalAboveThresh* diffLogic = new DifferentiateIfSignalAboveThresh(new_cell.get(), startingState);
         diffLogic->SetThresh(20);
         new_cell->SetLogic<Differentiation>(diffLogic);

         cells.push_back(new_cell);
      }

      // Use variable radii 
      NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(*mesh, cells);
      cellPopulation->SetUseVariableRadii(true);

      // Add boundary conditions forming a column, with a killer at the top.
      std::vector< boost::shared_ptr<PlaneBoundaryCondition<3> > > walls
                    = GetBoxWalls(cellPopulation, sideOfInitialCellCube);
      boost::shared_ptr<PlaneBasedCellKiller<3> > killer
                    = GetTopKiller(cellPopulation, sideOfInitialCellCube);

      OffLatticeSimulation<3>* simulation = new OffLatticeSimulation<3>(*cellPopulation);
      simulation->AddCellPopulationBoundaryCondition(walls[0]);
      simulation->AddCellPopulationBoundaryCondition(walls[1]);
      simulation->AddCellPopulationBoundaryCondition(walls[2]);
      simulation->AddCellPopulationBoundaryCondition(walls[3]);
      simulation->AddCellPopulationBoundaryCondition(walls[4]);
      simulation->AddCellKiller(killer);

      // Add a niche signal. Oriented along z axis, originates at z=0.
      MAKE_PTR(PositionSignal<3>, p_nicheSignal);
      p_nicheSignal->SetAxisOfInterest(2);
      p_nicheSignal->SetMeasurementStartPoint(0);
      simulation->AddSimulationModifier(p_nicheSignal);

      // Add a simple repulsion force
      MAKE_PTR(RepulsionForce<3>, p_force);
      simulation->AddForce(p_force);

      // Set up simulation and run
      simulation->SetOutputDirectory("TestLogicCellColumn");
      simulation->SetSamplingTimestepMultiple(10);
      simulation->SetDt(1.0/120.0);
      simulation->SetEndTime(30);
      simulation->Solve();
   }

};

#endif //TESTLOGICCELLPOPULATION_HPP