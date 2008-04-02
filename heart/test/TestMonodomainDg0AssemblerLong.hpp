/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _TESTMONODOMAINDG0ASSEMBLERLONG_HPP_
#define _TESTMONODOMAINDG0ASSEMBLERLONG_HPP_



#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petscvec.h>
#include <vector>

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ReplicatableVector.hpp"
#include "CheckMonoLr91Vars.hpp"


#include <time.h>

class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    InitialStimulus *mpStimulus;
    unsigned mNodeNum;
public:
    PointStimulus2dCellFactory(unsigned nodeNum)
            : AbstractCardiacCellFactory<2>(0.01),
            mNodeNum(nodeNum)
    {
        mpStimulus = new InitialStimulus(-6000.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node == mNodeNum)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        }
    }
    
    ~PointStimulus2dCellFactory(void)
    {
        delete mpStimulus;
    }
};


class TestMonodomainDg0AssemblerLong : public CxxTest::TestSuite
{
public:

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    // We run for 500 ms and then check that all the voltages at the final time
    // have returned to the resting potential of -84.5
    // test should take about 30mins (or less)
    void TestMonodomainDg02DWithPointStimulusInTheVeryCentreOfTheMesh( void )
    {
        
        PointStimulus2dCellFactory cell_factory(60); // Central node
        
        MonodomainProblem<2> monodomain_problem(&cell_factory);
        
        monodomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomain_problem.SetEndTime(500);   // 500 ms
        monodomain_problem.SetOutputDirectory("MonoDg02dWithPointStimulusLong");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainLR91_2dWithPointStimulusLong");
        monodomain_problem.SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        monodomain_problem.Initialise();
        
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);        
        
        monodomain_problem.Solve();
        
        CheckMonoLr91Vars(monodomain_problem);
        
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        
        /*
        * Test that corners are 'equal', and centres of sides.
        * Irregularities in which way the triangles are oriented make
        * this rather difficult, especially since the edges are sampled
        * during the upstroke.
        */
        
        // corners
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[10],  0.1);
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[110], 0.1);
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[120], 0.1);
        
        // centres of edges
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[55],  0.1);
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[65],  0.1);
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[115], 0.1);
        
        int num_nodes = monodomain_problem.rGetMesh().GetNumNodes();
        // test final voltages have returned to the resting potential
        for (int i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(voltage_replicated[i], -84.5, 1);
        }
        
    }
};
#endif //_TESTMONODOMAINDG0ASSEMBLERLONG_HPP_
