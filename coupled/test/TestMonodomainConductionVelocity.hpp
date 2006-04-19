#ifndef _TESTMONODOMAINCONDUCTIONVELOCITY_HPP_
#define _TESTMONODOMAINCONDUCTIONVELOCITY_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"


#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
//#include <iostream>

#include "ConformingTetrahedralMesh.cpp"
#include "PropagationPropertiesCalculator.hpp"
#include "ColumnDataReader.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblemIteration7.hpp"
#include "AbstractCardiacCellFactory.hpp"


class PointStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PointStimulusCellFactory(double timeStep) : AbstractCardiacCellFactory<1>(timeStep)
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(-600, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (node == 0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus, mTimeStep);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
        }
        
    }
    
    ~PointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestMonodomainConductionVelocity : public CxxTest::TestSuite 
{
public:
    // Solve on a 1D string of cells, 1cm long with a space step of 0.1mm.
    void TestMonodomainDg01D_100elements()
    {
        PointStimulusCellFactory cell_factory(0.01); // ODE time step (ms)
        MonodomainProblemIteration7<1> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        monodomain_problem.SetEndTime(30);   // 30 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg01d");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");

        monodomain_problem.Initialise();
        monodomain_problem.Solve();
        
        double* voltage_array;
    
        // test whether voltages and gating variables are in correct ranges

        VecGetArray(monodomain_problem.mCurrentVoltage, &voltage_array); 
        
        for(int global_index=monodomain_problem.mLo; global_index<monodomain_problem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomain_problem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomain_problem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.mMonodomainPde->GetCardiacCell(global_index)->GetStateVariables();
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }   
            }
        }
        VecRestoreArray(monodomain_problem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomain_problem.mCurrentVoltage);
        VecAssemblyEnd(monodomain_problem.mCurrentVoltage);
        VecDestroy(monodomain_problem.mCurrentVoltage);
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            // Calculate the conduction velocity
            ColumnDataReader simulation_data("testoutput/MonoDg01d",
                                             "NewMonodomainLR91_1d");
            PropagationPropertiesCalculator ppc(&simulation_data);
            double velocity;
            
            // Check action potential propagated to node 95
            TS_ASSERT_THROWS_NOTHING(velocity=ppc.CalculateConductionVelocity(5,95,0.9));
            
            // The value should be approximately 50cm/sec
            // i.e. 0.05 cm/msec (which is the units of the simulation)
            TS_ASSERT_DELTA(velocity, 0.05, 0.003);
        }
    }

    // Solve on a 1D string of cells, 1cm long with a space step of 0.05mm.
    // The purpose of this test is to check that it still gives similar results
    // to TestMonodomainDg01D_100elements. This is based on the assumption that 
    // 0.1mm is the space step from which convergence occurs. In that context, 
    // going for a smaller space step should still the same converging results. 
    // Note though that this requires a smaller time step for the integrator...
    void TestMonodomainDg01D_200elements()
    {
        PointStimulusCellFactory cell_factory(0.002);  // ODE time step (ms)
        MonodomainProblemIteration7<1> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_200_elements");
        monodomain_problem.SetEndTime(30);   // 30 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg01d");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");
        monodomain_problem.SetPdeTimeStep(0.002); // ms

        monodomain_problem.Initialise();
        monodomain_problem.Solve();
        
        double* voltage_array;
    
        // test whether voltages and gating variables are in correct ranges

        VecGetArray(monodomain_problem.mCurrentVoltage, &voltage_array); 
        
        for(int global_index=monodomain_problem.mLo; global_index<monodomain_problem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomain_problem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomain_problem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.mMonodomainPde->GetCardiacCell(global_index)->GetStateVariables();
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }   
            }
        }
        VecRestoreArray(monodomain_problem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomain_problem.mCurrentVoltage);
        VecAssemblyEnd(monodomain_problem.mCurrentVoltage);
        VecDestroy(monodomain_problem.mCurrentVoltage);
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            // Calculate the conduction velocity
            ColumnDataReader simulation_data("testoutput/MonoDg01d",
                                             "NewMonodomainLR91_1d");
            PropagationPropertiesCalculator ppc(&simulation_data);
            double velocity;
            
            // Check action potential propagated to node 185
            TS_ASSERT_THROWS_NOTHING(velocity=ppc.CalculateConductionVelocity(5,185,0.9));
            
            // The value should be approximately 50cm/sec
            // i.e. 0.05 cm/msec (which is the units of the simulation)
            TS_ASSERT_DELTA(velocity, 0.05, 0.003);
        }
    }

    // Solve on a 1D string of cells, 1cm long with a space step of 0.5mm.
    // See above for the reasons behind this test. 
    // Note that this space step ought to be too big!
    void TestMonodomainDg01D_20elements()
    {
        PointStimulusCellFactory cell_factory(0.01); // ODE time step (ms)
        MonodomainProblemIteration7<1> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_20_elements");
        monodomain_problem.SetEndTime(30);   // 30 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg01d");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");
        monodomain_problem.Initialise();
        monodomain_problem.Solve();
        
        double* voltage_array;
    
        // test whether voltages and gating variables are in correct ranges

        VecGetArray(monodomain_problem.mCurrentVoltage, &voltage_array); 
        
        for(int global_index=monodomain_problem.mLo; global_index<monodomain_problem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomain_problem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomain_problem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.mMonodomainPde->GetCardiacCell(global_index)->GetStateVariables();
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }   
            }
        }
        VecRestoreArray(monodomain_problem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomain_problem.mCurrentVoltage);
        VecAssemblyEnd(monodomain_problem.mCurrentVoltage);
        VecDestroy(monodomain_problem.mCurrentVoltage);
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            // Calculate the conduction velocity
            ColumnDataReader simulation_data("testoutput/MonoDg01d",
                                             "NewMonodomainLR91_1d");
            PropagationPropertiesCalculator ppc(&simulation_data);
            double velocity;
            
            // Check action potential propagated to node 47
            TS_ASSERT_THROWS_NOTHING(velocity=ppc.CalculateConductionVelocity(1,19,0.9));
            
            // The value should NOT be approximately 50cm/sec
            // i.e. 0.05 cm/msec (which is the units of the simulation)
            TS_ASSERT(fabs(velocity-0.05) > 0.003);
        }
    }
};
#endif //_TESTMONODOMAINCONDUCTIONVELOCITY_HPP_
