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

#ifndef TESTNUMERICS_HPP_
#define TESTNUMERICS_HPP_

#include <cxxtest/TestSuite.h>

//Misc
#include "CellBasedEventHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "FixedDivisionTimingsCellCycleModel.hpp"
#include "SimpleCellCentrePositionTracker.hpp"
#include "CommandLineArguments.hpp"
#include "VoronoiDataWriter.hpp"

//Populations
#include "AbstractOffLatticeCellPopulation.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

//Forces etc
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PopulationTestingForce.hpp"

//Methods
#include "ForwardEulerNumericalMethod.hpp"
#include "BackwardEulerNumericalMethod.hpp"
#include "AdamsMoultonNumericalMethod.hpp"
#include "RK4NumericalMethod.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
* Tests different numerical method and cell population combinations. 
* Tracks the position of cell centres at regular intervals. 
*/

class TestNumerics : public AbstractCellBasedTestSuite
{

private:

	void ResetForNewRep(int rep){
        
        RandomNumberGenerator::Instance()->Reseed(rep);
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellId::ResetMaxCellId();
        CellBasedEventHandler::Reset();
    };


    void SetFilename(std::string& filename, std::string population, std::string cycle, std::string method, bool isAdaptive, int power, int rep){

    	std::stringstream repAsString;
		repAsString << rep;
		std::stringstream powAsString;
		powAsString << "p" << power << "_";
		
		filename = population + "_" + cycle + "_" + method + "_";
		if(isAdaptive){
			filename += "true_";
		}else{
			filename += "false_"; 
		}
		filename += powAsString.str();
		filename += repAsString.str();
    };


    template<unsigned DIM>
    boost::shared_ptr<AbstractNumericalMethod<DIM> > MakeNumericalMethod(std::string method){
    	if(method=="FE"){
    		MAKE_PTR(ForwardEulerNumericalMethod<DIM>, nm);
    		return nm;
		}else if(method=="RK4"){
			MAKE_PTR(RK4NumericalMethod<DIM>, nm);
    		return nm;
		}else if(method=="BE"){
			MAKE_PTR(BackwardEulerNumericalMethod<DIM>, nm);
    		return nm;
		}else if(method=="AM2"){
			MAKE_PTR(AdamsMoultonNumericalMethod<DIM>, nm);
    		return nm;
		}else{
			EXCEPTION("Unrecognised numerical method");
		}
		MAKE_PTR(ForwardEulerNumericalMethod<DIM>, defaultNm);
		return defaultNm;
    }


    template<unsigned DIM>
	void GenerateStemCellsWithDistributedAges(std::vector<CellPtr>& rCells,
                                              const std::vector<unsigned> realCellIndices,
                                              std::string cycle,
                                              double ccLength)
	{
    	assert(!realCellIndices.empty());
    	unsigned nCells = realCellIndices.size();

    	rCells.clear();
    	rCells.reserve(nCells);
    	CellPropertyRegistry::Instance()->Clear();

    	for (unsigned i=0; i<nCells; i++)
    	{	
    		AbstractCellCycleModel* pCCM;
    		if(cycle=="stoch"){
    			pCCM = new StochasticDurationCellCycleModel();
    			pCCM->SetStemCellG1Duration(0.583*ccLength);
    			if(0.583*ccLength < 2){
    				EXCEPTION("G1 duration < 2 hours, the default range of stochastic variation.");
    			}
    			pCCM->SetG2Duration(0.167*ccLength);
    			pCCM->SetMDuration(0.042*ccLength);
    			pCCM->SetSDuration(0.208*ccLength);
    		}else{
    			pCCM = new FixedDivisionTimingsCellCycleModel(ccLength);
    		}
    		pCCM->SetDimension(DIM);

    	    boost::shared_ptr<AbstractCellProperty> pState(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
    	    CellPtr pCell(new Cell(pState, pCCM));

			pCell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

    	    double birthTime = - ccLength*RandomNumberGenerator::Instance()->ranf();
    	    pCell->SetBirthTime(birthTime);
    	    rCells.push_back(pCell);
    	}
	}


	template<unsigned DIM>
	void SetAMT(OffLatticeSimulation<DIM>* simulation, float AMT){
		AbstractCellPopulation<DIM,DIM>* pop = &(simulation->rGetCellPopulation());
        dynamic_cast<AbstractOffLatticeCellPopulation<DIM,DIM>*>(pop)->SetAbsoluteMovementThreshold((double)AMT);
	}


	template<unsigned DIM>
    void AddTracking(OffLatticeSimulation<DIM>* simulation, int interval){
    	MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<DIM>, trackingModifier,(interval,1));
    	simulation->AddSimulationModifier(trackingModifier);
    }


    template<unsigned DIM>
    void AddForce(OffLatticeSimulation<DIM>* simulation, std::string population, double cutoff){
    	
    	if(population=="node" || population=="mesh" || population=="ghost"){
    		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, pForce);
    		pForce->SetCutOffLength(cutoff);
    		simulation->AddForce(pForce);

    	}else if(population=="particle"){
    		MAKE_PTR(PopulationTestingForce<DIM>, pForce);
    		simulation->AddForce(pForce);

    	}else if(population=="vertex"){
    		MAKE_PTR(NagaiHondaForce<DIM>, pForce);
            simulation->AddForce(pForce);
            MAKE_PTR(SimpleTargetAreaModifier<DIM>, pGrowthModifier);
            simulation->AddSimulationModifier(pGrowthModifier);
    	}
    }


    void CleanupNodes(std::vector< Node<3>* > nodes)
    {
    	for(int i=0; i<nodes.size(); i++){
    		delete nodes[i];
    	}
    };


public:

	void TestWithPositionRecording() throw (Exception)
	{
		int* argcount = CommandLineArguments::Instance()->p_argc;
		TS_ASSERT(*argcount==10);
		
		char*** args = CommandLineArguments::Instance()->p_argv;
		std::string population = (*args)[1];
		std::string cycle = (*args)[2];	
		std::string method = (*args)[3];
		int adaptive = atoi((*args)[4]);	
		float AMT = atof((*args)[5]);
        int stepPowerMin = atoi((*args)[6]);
        int stepPowerMax = atoi((*args)[7]);
        int rangelower = atoi((*args)[8]);
        int rangeupper = atoi((*args)[9]);
		bool isAdaptive = (adaptive==1); 

		// SETUP CONSTANTS
		double endTime = 30;
		double interactionCutoff = 3.1;
		double ccLength = 6;
		int nCellsPerSide = 2;
    	int nCells = nCellsPerSide*nCellsPerSide;
        double dampingNormal = 1.1;

		for(int power = stepPowerMin; power <= stepPowerMax; power++){

			for(int rep = rangelower; rep <= rangeupper; rep++){

				ResetForNewRep(rep);
				std::string filename;
				SetFilename(filename, population, cycle, method, isAdaptive, power, rep);
				std::cout << filename << std::endl;

				double dt = 1.0/((double)pow(2,power));
				

				// =================================================================================================== 
				// 3D SIMULATION ===================================================================================== 
				// =================================================================================================== 
				if(population=="node" || population=="particle"){

					boost::shared_ptr<AbstractNumericalMethod<3> > numericalMethod = MakeNumericalMethod<3>(method);
					
					OffLatticeSimulation<3>* simulation;
					std::vector<Node<3>*> nodes;
        			std::vector<unsigned> realCellIndices;
        			std::vector<CellPtr> cells;
        			AbstractOffLatticeCellPopulation<3,3>* cellPopulation;

        			if(population=="node"){ 

        				for(int n=0; n<nCells; n++){
        					nodes.push_back(new Node<3>((unsigned)n, false, n/nCellsPerSide, n%nCellsPerSide, 0.0));
        					realCellIndices.push_back((unsigned)n);
        				}

        			}else{

        				for(int n=0; n<nCells; n++){
            				if(n<nCells/2){
            					nodes.push_back(new Node<3>((unsigned)n, false, n/nCellsPerSide, n%nCellsPerSide, 0.0));
            					realCellIndices.push_back((unsigned)n);
            				}else{
            					nodes.push_back(new Node<3>((unsigned)n, true, n/nCellsPerSide, n%nCellsPerSide, 0.0));
            				}
            			}
            		}

        			MAKE_PTR(NodesOnlyMesh<3>, pMesh);
        			pMesh->ConstructNodesWithoutMesh(nodes, interactionCutoff);
        			GenerateStemCellsWithDistributedAges<3>(cells, realCellIndices, cycle, ccLength);       
        			
        			if(population=="node"){
        				cellPopulation = new NodeBasedCellPopulation<3>(*pMesh, cells);
        			}else{
        				cellPopulation = new NodeBasedCellPopulationWithParticles<3>(*pMesh, cells, realCellIndices);
        			}

                    cellPopulation->SetDampingConstantNormal(dampingNormal);
        			simulation = new OffLatticeSimulation<3>(*cellPopulation, false, true, numericalMethod, isAdaptive);
					simulation->SetOutputDirectory(filename.c_str());
        	    	simulation->SetDt(dt);
        	    	simulation->SetSamplingTimestepMultiple((int)pow(2,power));
        	    	simulation->SetEndTime(endTime);

        	    	SetAMT<3>(simulation, AMT);
        	    	AddForce<3>(simulation, population, interactionCutoff);
        	    	AddTracking<3>(simulation, (int)pow(2,power));

        	    	simulation->Solve();

        	    	CleanupNodes(nodes);
        	    	delete simulation;
        	    	delete cellPopulation;
        	    



        	    // =================================================================================================== 
				// 2D SIMULATION ===================================================================================== 
				// =================================================================================================== 
				}else{

					boost::shared_ptr<AbstractNumericalMethod<2> > numericalMethod = MakeNumericalMethod<2>(method);
				    
				    OffLatticeSimulation<2>* simulation;
				    std::vector<unsigned> realCellIndices;
				    std::vector<CellPtr> cells;
				    HoneycombMeshGenerator* gen;
				    HoneycombVertexMeshGenerator* vGen;
				    AbstractOffLatticeCellPopulation<2,2>* cellPopulation;

				   	if(population=="mesh"){
		    			gen = new HoneycombMeshGenerator(nCellsPerSide, nCellsPerSide, 0);
		    			realCellIndices = gen->GetCellLocationIndices();
		    			GenerateStemCellsWithDistributedAges<2>(cells, realCellIndices, cycle, ccLength);  
            			cellPopulation = new MeshBasedCellPopulation<2>(*(gen->GetMesh()), cells);
            			dynamic_cast<MeshBasedCellPopulation<2>*>(cellPopulation)->CreateVoronoiTessellation();
                        dynamic_cast<MeshBasedCellPopulation<2>*>(cellPopulation)->AddPopulationWriter<VoronoiDataWriter>();
		    		
                    }else if(population=="ghost"){
		    			gen = new HoneycombMeshGenerator(nCellsPerSide, nCellsPerSide, 2);
		    			realCellIndices = gen->GetCellLocationIndices();
		    			GenerateStemCellsWithDistributedAges<2>(cells, realCellIndices, cycle, ccLength);  
  					    cellPopulation = new MeshBasedCellPopulationWithGhostNodes<2>(*(gen->GetMesh()), cells, realCellIndices);
  					    dynamic_cast<MeshBasedCellPopulation<2>*>(cellPopulation)->CreateVoronoiTessellation();
                        dynamic_cast<MeshBasedCellPopulation<2>*>(cellPopulation)->AddPopulationWriter<VoronoiDataWriter>();
		    		
                    }else{
		    			vGen = new HoneycombVertexMeshGenerator(nCellsPerSide, nCellsPerSide, false, 0.1);
		    			for(int i=0; i<nCells; i++){
  					    	realCellIndices.push_back((unsigned)i);
  					    }
  					    GenerateStemCellsWithDistributedAges<2>(cells, realCellIndices, cycle, ccLength);  
  					    cellPopulation = new VertexBasedCellPopulation<2>(*(vGen->GetMesh()), cells);
		    		}
     
                    cellPopulation->SetDampingConstantNormal(dampingNormal);
                    simulation = new OffLatticeSimulation<2>(*cellPopulation, false, true, numericalMethod, isAdaptive);
					simulation->SetOutputDirectory(filename.c_str());
        	    	simulation->SetDt(dt);
        	    	simulation->SetSamplingTimestepMultiple((int)pow(2,power));
        	    	simulation->SetEndTime(endTime);

        	    	SetAMT<2>(simulation, AMT);
        	    	AddForce<2>(simulation, population, interactionCutoff);
        	    	AddTracking<2>(simulation, (int)pow(2,power));

        	    	simulation->Solve();

        	    	delete simulation; 
        	    	delete cellPopulation;
        	    	if(population=="vertex"){
        	    		delete vGen;
        	    	}else{
        	    		delete gen;
        	    	}
				}
			}

		}
	};

};

#endif /*TESTNUMERICS_HPP_*/