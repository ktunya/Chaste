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

#include "LogicCell.hpp"
#include "EnvironmentalSignalTypes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellCycleWithDivisionMechanism.hpp"
#include "OffLatticeSimulation.hpp"
#include "RepulsionForce.hpp"
#include "SimpleAssymmetricStemDivision.hpp"
#include "DecayingNicheSignal.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "TempStateOutputViaCellData.hpp"
#include "DummyMutationState.hpp"
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

   void TestNodeBasedPopulationOfLogicCells(){

      //SETUP INITIAL CELL LOCATIONS
      int sideOfInitialCellCube = 3;
      int nCells = sideOfInitialCellCube * sideOfInitialCellCube * sideOfInitialCellCube;
      NodesOnlyMesh<3>* mesh  = make3DCubeMesh(sideOfInitialCellCube);

      //CREATE CELLS
      std::vector<CellPtr> cells;
      cells.clear();
      cells.reserve(nCells);

      for (int i=0; i < nCells; i++)
      {

         //MAKE A CELL
         MAKE_PTR(DummyMutationState, pDummyMut);
         MAKE_PTR_ARGS(LogicCell, new_cell, (pDummyMut, false));

         //ADD A CELL CYCLE
         CellCycleWithDivisionMechanism* cycleLogic = new CellCycleWithDivisionMechanism(new_cell.get(), CellCycle::G1, 0.0, 2.0);
         new_cell->setLogic<CellCycle>(cycleLogic);

         //ADD CELL GROWTH LOGIC
         CellGrowth* growthLogic = new CellGrowth(new_cell.get(), Growth::GROWING, 1.0, 0.5, 2.0);
         new_cell->setLogic<Growth>(growthLogic);

         //ADD DIFFERENTIATION LOGIC; START WITH SOME STEM AND SOME TA CELLS
         SimpleAssymmetricStemDivision* diffLogic;
         if(i < sideOfInitialCellCube * sideOfInitialCellCube){
            diffLogic = new SimpleAssymmetricStemDivision(new_cell.get(), Differentiation::S, 0);
         }else{
            diffLogic = new SimpleAssymmetricStemDivision(new_cell.get(), Differentiation::TA, 0);
         }
         new_cell->setLogic<Differentiation>(diffLogic);

         //DONE
         cells.push_back(new_cell);
      }

      //USE VARIABLE RADII
      NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(*mesh, cells);
      cellPopulation->SetUseVariableRadii(true);

      //SET BOUNDARY CONDITIONS AND TOP KILLER
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

      //ADD AN ENVIRONMENTAL NICHE SIGNAL
      MAKE_PTR(DecayingNicheSignal<3>, p_nicheSignal);
      simulation->AddSimulationModifier(p_nicheSignal);

      //ADD FORCE
      MAKE_PTR(RepulsionForce<3>, p_force);
      simulation->AddForce(p_force);

      //FINAL SIMULATION PROPERTIES AND SOLVE
      simulation->SetOutputDirectory("TestNodeBasedLogicCellPopulation");
      simulation->SetSamplingTimestepMultiple(10);
      simulation->SetDt(1.0/250.0);
      simulation->SetEndTime(30);
      simulation->Solve();
   }




};


#endif //TESTLOGICCELLPOPULATION_HPP
