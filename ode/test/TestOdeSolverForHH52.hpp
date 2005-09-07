#ifndef _TESTODESOLVERFORHH52_HPP_
#define _TESTODESOLVERFORHH52_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"
#include "OdeSolution.hpp"
#include "ColumnDataWriter.hpp"
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"


class TestOdeSolverForHH52 : public CxxTest::TestSuite
{
    public:
    
    
    // Test Ode Solver for HH52
    void testOdeSolverForHH52WithInitialStimulus(void)
    {
        /*
         * Set initial conditions and magnitude of stimulus
         * 
         */
		double voltage = 0.0;
        double n = 0.31768;
        double h = 0.59612;
        double m = 0.05293;
  
        double magnitudeOfStimulus = -20.0;  
        double durationOfStimulus  = 0.5 ;  // ms                     
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initialConditions;
        initialConditions.push_back(voltage);
        initialConditions.push_back(n);
        initialConditions.push_back(h);
        initialConditions.push_back(m);        
        /*
         * Choose function for stimulus
         */             
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus);
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(&stimulus);
        
        /*
         * Choose an ode solver
         */
        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1.0;
        double timeStep = 0.01;             
                
        OdeSolution solution_new = solver.Solve(&hh52_ode_system, startTime, endTime, timeStep, initialConditions);
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
		        
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("testoutput","HH52Result");
        int time_var_id=mpNewTestWriter->DefineFixedDimension("Time","ms", solution_new.mSolutions.size());
        int v_var_id = mpNewTestWriter->DefineVariable("V","milliamperes");
        int n_var_id = mpNewTestWriter->DefineVariable("n"," ");
        int h_var_id = mpNewTestWriter->DefineVariable("h"," ");
        int m_var_id = mpNewTestWriter->DefineVariable("m"," ");
        mpNewTestWriter->EndDefineMode();
				
        for (int i = 0; i < solution_new.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(time_var_id, solution_new.mTime[i], i);
            mpNewTestWriter->PutVariable(v_var_id, solution_new.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(n_var_id, solution_new.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(h_var_id, solution_new.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(m_var_id, solution_new.mSolutions[i][3], i);         
        }
        
        delete mpNewTestWriter;
        
        
        //read in good data file and compare line by line
        std::ifstream testfile("testoutput/HH52Result.dat",std::ios::in);
        std::ifstream goodfile("ode/test/data/HH52ResultGood.dat",std::ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(testfile, teststring))
        {
              getline(goodfile,goodstring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
                               
    }	
    
    void testOdeSolverForHH52WithRegularStimulus(void)
    {
        
        /*
         * This test doesn't really do anything. 
         * This is because it would have to run for quite a while and use up 
         * lots of memory to generate repeated stimulus results. To get lovely 
         * pictures increase the amount of time this runs for...
         * 
         */
        
        /*
         * Set initial conditions and magnitude of stimulus
         * 
         */
        double voltage = 0.0;
        double n = 0.31768;
        double h = 0.59612;
        double m = 0.05293;
  
        double magnitudeOfStimulus = -20.0;  
        double durationOfStimulus  = 0.5 ;  // ms                     
        double frequency = 1.0/50.0;
        double startStimulus = 40.0;              
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initialConditions;
        initialConditions.push_back(voltage);
        initialConditions.push_back(n);
        initialConditions.push_back(h);
        initialConditions.push_back(m);        
        /*
         * Choose function for stimulus
         */             
       
        RegularStimulus stimulus(magnitudeOfStimulus, durationOfStimulus, frequency, startStimulus);
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(&stimulus);
        
        /*
         * Choose an ode solver
         */      
        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 1.0;
        double timeStep = 0.01;             
                
        OdeSolution solution_new = solver.Solve(&hh52_ode_system, startTime, endTime, timeStep, initialConditions);
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
                
        ColumnDataWriter *mpNewTestWriter;
        mpNewTestWriter = new ColumnDataWriter("testoutput","HH52RegResult");
        int time_var_id=mpNewTestWriter->DefineFixedDimension("Time","ms", solution_new.mSolutions.size());
        int v_var_id = mpNewTestWriter->DefineVariable("V","milliamperes");
        int n_var_id = mpNewTestWriter->DefineVariable("n"," ");
        int h_var_id = mpNewTestWriter->DefineVariable("h"," ");
        int m_var_id = mpNewTestWriter->DefineVariable("m"," ");
        mpNewTestWriter->EndDefineMode();
                
        for (int i = 0; i < solution_new.mSolutions.size(); i++) 
        {
            mpNewTestWriter->PutVariable(time_var_id, solution_new.mTime[i], i);
            mpNewTestWriter->PutVariable(v_var_id, solution_new.mSolutions[i][0], i);
            mpNewTestWriter->PutVariable(n_var_id, solution_new.mSolutions[i][1], i);
            mpNewTestWriter->PutVariable(h_var_id, solution_new.mSolutions[i][2], i);
            mpNewTestWriter->PutVariable(m_var_id, solution_new.mSolutions[i][3], i);         
        }
        
        delete mpNewTestWriter;
        
        
//        //read in good data file and compare line by line
//        std::ifstream testfile("testoutput/NewLR91.dat",std::ios::in);
//        std::ifstream goodfile("ode/test/data/Lr91Good.dat",std::ios::in);
//        std::string teststring;
//        std::string goodstring;
//        while(getline(testfile, teststring))
//        {
//              getline(goodfile,goodstring);
//              TS_ASSERT_EQUALS(teststring,goodstring);
//        }
//        testfile.close();
//        goodfile.close();
                               
    }   

};



#endif //_TESTODESOLVERFORHH52_HPP_
