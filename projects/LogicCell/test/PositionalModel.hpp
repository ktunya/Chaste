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

//Generic Chaste headers
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "RandomNumberGenerator.hpp"
#include <limits>

//Extra mechanics classes
#include "DistalArmBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "RepulsionForceSizeCorrected.hpp"
#include "AppliedResistance.hpp"

//Data output
#include "DetailedCellTracker.hpp"
#include "CellCyclePhaseOutput.hpp"
#include "CellMitosisOutput.hpp"
#include "ZoneOutputOutput.hpp"

//Machinery for logic modules 
#include "LogicCell.hpp"
#include "AbstractCellLogic.hpp"
#include "LogicTypes.hpp"
#include "DummyMutationState.hpp"

//Our specific cell logic models, for growth, the cell cycle, and differentiation
#include "CellGrowthWithDiffDependence.hpp"
#include "StochasticGermCellCycleWithDiffGating.hpp"
#include "DifferentiateIfSignalAboveThresh.hpp"
#include "PositionSignal.hpp"

#ifndef POSITIONALMODEL_HPP
#define POSITIONALMODEL_HPP


class TestPositionalModel : public AbstractCellBasedTestSuite {

public:

	void TestPositional(){

        //Parameters------------------------------------------------------------------------------------------
        std::string outputDirectory( (*(CommandLineArguments::Instance()->p_argv))[1] ); 
        double resist =  atof((*(CommandLineArguments::Instance()->p_argv))[2]); 
        double spring =   atof((*(CommandLineArguments::Instance()->p_argv))[3]); 
        double zonelen = atof((*(CommandLineArguments::Instance()->p_argv))[4]); 
        double endtime = atof((*(CommandLineArguments::Instance()->p_argv))[5]); 
        int jobNumber = atoi((*(CommandLineArguments::Instance()->p_argv))[6]); 
        int stepsPerHour = atoi((*(CommandLineArguments::Instance()->p_argv))[7]); 

        std::cout << "JOBNUMBER: " << jobNumber << std::endl;

        int seed = (int)std::time(NULL);
        RandomNumberGenerator::Instance()->Reseed(seed); 

        double stoch = 0.2; 

		double length = 248;
		double radius = 11.3;
        double initialCellRowSpacing = 5;
        int    cellsPerRow = 10;

        //In hours; (G1, S, G2, M, MeioticG1, Meiotic S)
        int total = 8;
        double phases[] = {0.02*total, 0.57*total, 0.39*total,0.02*total,0.02*total,1.14*total};
        std::vector<double> cellCyclePhaseDurations(phases,  phases + sizeof(phases)/sizeof(double)); 

        //-----------------------------------------------------------------------------------------------------


		//Setup cells------------------------------------------------------------------------------------------
		
        NodesOnlyMesh<3>* mesh = CreateDistalArmCellLocations(length, radius, initialCellRowSpacing, cellsPerRow);
		
		//Create the cell population
		int nCells = mesh->GetNumNodes();
        int nProlifCells = cellsPerRow * (int)( zonelen / initialCellRowSpacing);

        std::vector<CellPtr> cells;
        cells.clear();
        cells.reserve(nCells);

        for (int i=0; i < nCells; i++)
        {
        	//Make a new empty logic cell
        	MAKE_PTR(DummyMutationState, pDummyMut);
        	MAKE_PTR_ARGS(LogicCell, new_cell, (pDummyMut, false));
            bool initiallyProliferative = (i < nProlifCells);

        	//Add a cell cycle with a random starting phase and a random elapsed time in phase. If the cell
            //is outside the proliferative region, start it in G0 with timeElapsed=doubleMax (i.e. no longer cycling).
            std::vector<double> initialCycleState;
            if(initiallyProliferative){
                initialCycleState = GetRandomCellCycleStart(cellCyclePhaseDurations);
            }else{
                initialCycleState.push_back(CellCycle::G0);
                initialCycleState.push_back(std::numeric_limits<double>::max());
                new_cell->GetCellData()->SetItem("TimeAtHim3Labelling",0);
            }

            StochasticGermCellCycleWithDiffGating* cycleLogic = new StochasticGermCellCycleWithDiffGating(
                                                                    new_cell.get(), (int)initialCycleState[0], initialCycleState[1]);
            cycleLogic->setPhaseDurations(cellCyclePhaseDurations);
            cycleLogic->setStochasticity(stoch);
            cycleLogic->setGate(CellCycle::G1);
            new_cell->setLogic<CellCycle>(cycleLogic);


            //Add a differentiation logic. Start cells in the stem state within the proliferative region 
            //and the diff state outside it
            int initialDiffState;
            if(initiallyProliferative){
                initialDiffState = Differentiation::S;
            }else{
                initialDiffState = Differentiation::D;
            }
            DifferentiateIfSignalAboveThresh* diffLogic = new DifferentiateIfSignalAboveThresh(new_cell.get(), initialDiffState);
            diffLogic->setThreshSignal(zonelen); 
            new_cell->setLogic<Differentiation>(diffLogic);
   

            //Add a cell growth logic
            CellGrowthWithDiffDependence* growthLogic = new CellGrowthWithDiffDependence(new_cell.get(), Growth::GROWING, 2.8);
            growthLogic->setHourlyRadialGrowthRate(0.04);
            growthLogic->setMaxRadiusDiffd(4.0);
            growthLogic->setMaxRadiusUndiffd(2.8);
            new_cell->setLogic<Growth>(growthLogic);
   
            cells.push_back(new_cell);
        }
        NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(*mesh, cells);
        cellPopulation->SetUseVariableRadii(true);
        cellPopulation->SetAbsoluteMovementThreshold(2.0);
        cellPopulation->SetDampingConstantNormal(1.0);
        //---------------------------------------------------------------------------------------------------


        //Setup rest of system-------------------------------------------------------------------------------

        //Apply the boundary condition
        MAKE_PTR_ARGS(DistalArmBoundaryCondition<3>, p_boundary, (cellPopulation, radius, length));
        OffLatticeSimulation<3>* simulation = new OffLatticeSimulation<3>(*cellPopulation);
        simulation->AddCellPopulationBoundaryCondition(p_boundary);

        //Add a plane based killer at the end of the arm
        c_vector<double, 3> point;
        point[0] = length; point[1] = 0; point[2] = 0;
        c_vector<double, 3> normal;
        normal[0] = 1; normal[1] = 0; normal[2] = 0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_killer, (cellPopulation, point, normal));
        simulation->AddCellKiller(p_killer);

        //Add a resistence force at end of distal arm
        double strength = resist;
        c_vector<double,3> direction;
        direction[0] = -1; direction[1] = 0; direction[2] = 0; 
        c_vector<double,3> targetRegionBottom;
        targetRegionBottom[0] = length-2.1; 
        targetRegionBottom[1] = -radius; 
        targetRegionBottom[2] = -radius;
        c_vector<double,3> targetRegionTop;
        targetRegionTop[0] = length; 
        targetRegionTop[1] = radius; 
        targetRegionTop[2] = radius;

        MAKE_PTR_ARGS(AppliedResistance<3>, p_resist,(strength, targetRegionBottom, targetRegionTop, direction));
        simulation->AddForce(p_resist);

        //Add cell-cell repulsion force
        MAKE_PTR(RepulsionForceSizeCorrected<3>, p_force);
        p_force->SetMeinekeSpringStiffness(spring);
        p_force->SetMeinekeSpringGrowthDuration(cellCyclePhaseDurations[3]);
        simulation->AddForce(p_force);

        //Add a positional signal from the DTC
        MAKE_PTR(PositionSignal<3>, p_signal);
        p_signal->SetAxisOfInterest(0);
        p_signal->SetMeasurementStartPoint(-radius);
        simulation->AddSimulationModifier(p_signal);

        //--------------------------------------------------------------------------------------------------

        //Add simulation output
        MAKE_PTR_ARGS(DetailedCellTracker<3>, p_detailedRecording, (stepsPerHour/2,5));
        simulation->AddSimulationModifier(p_detailedRecording);
        MAKE_PTR_ARGS(CellMitosisOutput<3>, p_mitosisRecording, (stepsPerHour/10, -radius, length, initialCellRowSpacing));
        simulation->AddSimulationModifier(p_mitosisRecording);
        MAKE_PTR_ARGS(CellCyclePhaseOutput<3>, p_phaseRecording, (stepsPerHour/10));
        simulation->AddSimulationModifier(p_phaseRecording);
        MAKE_PTR_ARGS(ZoneOutputOutput<3>, p_outputRecording, (stepsPerHour));
        simulation->AddSimulationModifier(p_outputRecording);

        //Set simulation properties and solve
        simulation->SetOutputDirectory(outputDirectory.c_str());
        simulation->SetSamplingTimestepMultiple(stepsPerHour);
        simulation->SetDt(1.0/(double)stepsPerHour);
        simulation->SetEndTime(endtime);
        simulation->Solve();

	}; 






private:

    NodesOnlyMesh<3>* CreateDistalArmCellLocations(double armLength, double armRadius, double rowSpacing, int cellsPerRow){

        double angleSpacing = 2*M_PI/cellsPerRow;
        
        std::vector< Node<3>* > nodes;
        unsigned nodeIndex = 0;

        for(double x = -armRadius; x < armLength; x += rowSpacing){
            for(double angle = 0; angle < 2*M_PI; angle += angleSpacing){
                Node<3>* n = new Node<3>(nodeIndex, false, x, armRadius*(cos(angle)), armRadius*(sin(angle)));
                n->SetRadius(2.8);
                nodes.push_back(n);
                nodeIndex++;
            }
        }

        NodesOnlyMesh<3>* distalMesh = new NodesOnlyMesh<3>;
        distalMesh->ConstructNodesWithoutMesh(nodes, 10);
        for (unsigned i = 0; i<nodes.size(); i++)
        {
           delete nodes[i];
        }
        return distalMesh;
    }


    std::vector<double> GetRandomCellCycleStart(std::vector<double> durations){

        int state;
        double timeInState = 0;

        double totalCycleLength = 0;
        for(int i=0; i<4; i++){
            totalCycleLength += durations[i];
        }

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance(); 
        double random = totalCycleLength*(p_gen->ranf());

        if(random < durations[0]){
            state = 0;
            timeInState = random;
        }else if(random >= durations[0] && random < (durations[0]+durations[1])){
            state = 1;
            timeInState = random - durations[0];
        }else if(random >= (durations[0]+durations[1]) && random < (durations[0]+durations[1]+durations[2])){
            state = 2;
            timeInState = random - durations[0]-durations[1];
        }else if(random >= (durations[0]+durations[1]+durations[2]) ){
            state = 3;
            timeInState = random - durations[0]-durations[1]-durations[2];
        }

        std::vector<double> result;
        result.push_back((double)state);
        result.push_back(timeInState);
        return(result);
    }
};

#endif //TESTPOSITIONAL_HPP