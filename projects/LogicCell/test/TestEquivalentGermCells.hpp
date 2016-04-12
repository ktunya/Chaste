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

// General Chaste headers
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "PlaneBasedCellKiller.hpp"
#include <limits>

// Extra mechanics simulation classes
#include "DistalArmBoundaryCondition.hpp"
#include "RepulsionForceSizeCorrected.hpp"
#include "AppliedResistance.hpp"

// Required for using logic modules 
#include "LogicCell.hpp"
#include "AbstractCellLogic.hpp"
#include "LogicTypes.hpp"
#include "DummyMutationState.hpp"

// Our specific intracellular logic for growth, the cell cycle, and differentiation
#include "CellGrowthWithDiffDependence.hpp"
#include "StochasticGermCellCycleWithDiffGating.hpp"
#include "DifferentiateIfSignalAboveThresh.hpp"
#include "PositionSignal.hpp"

// Data output
#include "DetailedCellTracker.hpp"
#include "CellCyclePhaseOutput.hpp"
#include "CellMitosisOutput.hpp"
#include "ZoneOutputOutput.hpp"

#ifndef POSITIONALMODEL_HPP
#define POSITIONALMODEL_HPP


class TestPositionalModel : public AbstractCellBasedTestSuite {

private: 

    void ResetForNewRun(){
        int seed = (int)std::time(NULL);
        RandomNumberGenerator::Instance()->Reseed(seed); 
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellId::ResetMaxCellId();
        CellBasedEventHandler::Reset();
    }


    NodesOnlyMesh<3>* CreateDistalArmCellLocations(double armLength, 
                                                   double armRadius, 
                                                   double rowSpacing, 
                                                   int cellsPerRow)
    {
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
        distalMesh->ConstructNodesWithoutMesh(nodes, 9);
        for (unsigned i = 0; i<nodes.size(); i++)
        {
           delete nodes[i];
        }
        return distalMesh;
    }


    std::vector<double> GetRandomCellCycleStart(std::vector<double> phaseDurations){

        int    state;
        double phaseStartTime = 0;

        double totalMitoticCellCycleLength = 0;
        for(int i=0; i<4; i++){
            totalMitoticCellCycleLength += phaseDurations[i];
        }

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance(); 
        double random = totalMitoticCellCycleLength*(p_gen->ranf());

        if(random < phaseDurations[0]){
            state = 0;
            phaseStartTime = -random;
        }else if(random >= durations[0] && random < (durations[0]+durations[1])){
            state = 1;
            phaseStartTime = -random + durations[0];
        }else if(random >= (durations[0]+durations[1]) && random < (durations[0]+durations[1]+durations[2])){
            state = 2;
            phaseStartTime = -random + durations[0]+durations[1];
        }else if(random >= (durations[0]+durations[1]+durations[2]) ){
            state = 3;
            phaseStartTime = -random + durations[0]+durations[1]+durations[2];
        }

        std::vector<double> result;
        result.push_back((double)state);
        result.push_back(phaseStartTime);
        return(result);
    }

    void SetFilename(std::string& filename, double springStrength, double resistanceForce, double zoneLengthX, int rep){

        std::stringstream fileNameStream;
        fileNameStream << springStrength << "_" << resistanceForce << "_" << zoneLengthX << "_" << rep;
        filename += fileNameStream.str();
    };

public:

	void TestEquivalentGermCellsModel(){

        //Parameters---------------------------------------------------------------------------------- 
        double resistanceForce =  atof((*(CommandLineArguments::Instance()->p_argv))[1]); 
        double springStrength =   atof((*(CommandLineArguments::Instance()->p_argv))[2]); 
        double zoneLengthX = atof((*(CommandLineArguments::Instance()->p_argv))[3]); 
        double endTime = atof((*(CommandLineArguments::Instance()->p_argv))[4]); 
        int stepsPerHour = atoi((*(CommandLineArguments::Instance()->p_argv))[5]); 
        int repStart = atoi((*(CommandLineArguments::Instance()->p_argv))[6]); 
        int repEnd = atoi((*(CommandLineArguments::Instance()->p_argv))[7]); 
        
        double stoch = 0.1; 
		double armLength = 248;
		double armRadius = 11.3;
        double initialCellRowSpacing = 5;
        int    initialCellsPerRow = 10;

        //Phase lengths in hours (G1, S, G2, M, MeioticG1, Meiotic S)
        int total = 8;
        std::vector<double> cellCyclePhaseDurations;
        cellCyclePhaseDurations.push_back(0.02*total);
        cellCyclePhaseDurations.push_back(0.57*total);
        cellCyclePhaseDurations.push_back(0.39*total);
        cellCyclePhaseDurations.push_back(0.02*total);
        cellCyclePhaseDurations.push_back(0.02*total);
        cellCyclePhaseDurations.push_back(1.14*total);
        //--------------------------------------------------------------------------------------------


        for(int rep = repStart; rep <= repEnd; rep++){

            ResetForNewRun();

		    //Set up cells--------------------------------------------------------------------------------
            NodesOnlyMesh<3>* mesh = CreateDistalArmCellLocations(armLength, 
                                                                  armRadius, 
                                                                  initialCellRowSpacing, 
                                                                  cellsPerRow);
		  
		    int nCells = mesh->GetNumNodes();
            int nProlifCells = initialCellsPerRow * (int)((zoneLengthX-armRadius)/initialCellRowSpacing);

            std::vector<CellPtr> cells;
            cells.clear();
            cells.reserve(nCells);

            for (int i=0; i<nCells; i++)
            {
            	// Make an empty logic cell
            	MAKE_PTR(DummyMutationState, pDummyMut);
            	MAKE_PTR_ARGS(LogicCell, new_cell, (pDummyMut, false));

                bool initiallyProliferative = (i<nProlifCells);
                
            	// Add a cell cycle with a random starting phase. If the cell is not initially proliferative
                // start it in G0 with phaseStartTime = -doubleMax.
                std::vector<double> startingCellCycleState;
                if(initiallyProliferative){
                    startingCellCycleState = GetRandomCellCycleStart(cellCyclePhaseDurations);
                }else{
                    startingCellCycleState.push_back(CellCycle::G0);
                    startingCellCycleState.push_back(std::numeric_limits<double>::max());
                    new_cell->GetCellData()->SetItem("TimeAtHim3Labelling",0);
                }

                StochasticGermCellCycleWithG1Diff* cycleLogic = new StochasticGermCellCycleWithG1Diff(new_cell.get(), 
                                                                                                      (int)startingCellCycleState[0],
                                                                                                      startingCellCycleState[1]);
                cycleLogic->SetPhaseDurations(cellCyclePhaseDurations);
                cycleLogic->SetStochasticity(stoch);
                new_cell->SetLogic<CellCycle>(cycleLogic);

                // Add differentiation logic. Cells in the proliferative region are S by virtue of position, 
                // while outside the zone become D
                int initialDiffState;
                if(initiallyProliferative){
                    initialDiffState = Differentiation::S;
                }else{
                    initialDiffState = Differentiation::D;
                }
                DifferentiateIfSignalAboveThresh* diffLogic = new DifferentiateIfSignalAboveThresh(new_cell.get(), initialDiffState);
                diffLogic->setThresh(zoneLengthX); 
                new_cell->SetLogic<Differentiation>(diffLogic);
   
                // Add linear cell growth, with a larger max radius for differentiated cells
                CellGrowthWithDiffDependence* growthLogic = new CellGrowthWithDiffDependence(new_cell.get(), Growth::GROWING, 2.8);
                growthLogic->setHourlyRadialGrowthRate(0.1);
                growthLogic->setMaxRadiusDiffd(4.0);
                growthLogic->setMaxRadiusUndiffd(2.8);
                new_cell->SetLogic<Growth>(growthLogic);
   
                cells.push_back(new_cell);
            }
            NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(*mesh, cells);
            cellPopulation->SetUseVariableRadii(true);
            cellPopulation->SetAbsoluteMovementThreshold(2.0);
            cellPopulation->SetDampingConstantNormal(1.0);
            //-------------------------------------------------------------------------------------------


            //Set up rest of system----------------------------------------------------------------------

            // Apply the boundary condition
            MAKE_PTR_ARGS(DistalArmBoundaryCondition<3>, p_boundary, (cellPopulation, armRadius, armLength));
            OffLatticeSimulation<3>* simulation = new OffLatticeSimulation<3>(*cellPopulation);
            simulation->AddCellPopulationBoundaryCondition(p_boundary);

            // Add a plane based cell killer at the end of the distal arm
            c_vector<double, 3> point;
            point[0] = armLength; point[1] = 0; point[2] = 0;
            c_vector<double, 3> normal;
            normal[0] = 1; normal[1] = 0; normal[2] = 0;
            MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_killer, (cellPopulation, point, normal));
            simulation->AddCellKiller(p_killer);

            // Add a resistance force at end of distal arm
            c_vector<double,3> direction;
            direction[0] = -1; direction[1] = 0; direction[2] = 0; 
            c_vector<double,3> targetRegionBottom;
            targetRegionBottom[0] = armLength-5.6; 
            targetRegionBottom[1] = -armRadius; 
            targetRegionBottom[2] = -armRadius;
            c_vector<double,3> targetRegionTop;
            targetRegionTop[0] = armLength; 
            targetRegionTop[1] = armRadius; 
            targetRegionTop[2] = armRadius;
            MAKE_PTR_ARGS(AppliedResistance<3>, p_resist,(resistanceForce, targetRegionBottom, targetRegionTop, direction));
            simulation->AddForce(p_resist);

            // Add cell-cell repulsion force
            MAKE_PTR(RepulsionForceSizeCorrected<3>, p_force);
            p_force->SetMeinekeSpringStiffness(springStrength);
            p_force->SetMeinekeSpringGrowthDuration(cellCyclePhaseDurations[3]);
            simulation->AddForce(p_force);

            // Add a position signal from the DTC
            MAKE_PTR(PositionSignal<3>, p_signal);
            p_signal->SetAxisOfInterest(0);
            p_signal->SetMeasurementStartPoint(-armRadius);
            simulation->AddSimulationModifier(p_signal);

            //--------------------------------------------------------------------------------------

            //Add simulation output
            MAKE_PTR_ARGS(DetailedCellTracker<3>, p_detailedRecording, (stepsPerHour/2,5));
            simulation->AddSimulationModifier(p_detailedRecording);
            MAKE_PTR_ARGS(CellMitosisOutput<3>, p_mitosisRecording, (stepsPerHour/10, -radius, length, initialCellRowSpacing));
            simulation->AddSimulationModifier(p_mitosisRecording);
            MAKE_PTR_ARGS(CellCyclePhaseOutput<3>, p_phaseRecording, (stepsPerHour/10));
            simulation->AddSimulationModifier(p_phaseRecording);
            MAKE_PTR_ARGS(ZoneOutputOutput<3>, p_outputRecording, (stepsPerHour));
            simulation->AddSimulationModifier(p_outputRecording);

            // Set simulation properties and solve
            std::string outputDirectory("EquivalentGermCells");
            SetFilename(outputDirectory, springStrength, resistanceForce, zoneLengthX, rep);
            simulation->SetOutputDirectory(outputDirectory.c_str());
            simulation->SetSamplingTimestepMultiple(stepsPerHour);
            simulation->SetDt(1.0/(double)stepsPerHour);
            simulation->SetEndTime(endTime);
            simulation->Solve();
        };
	};
};

#endif //TESTPOSITIONAL_HPP