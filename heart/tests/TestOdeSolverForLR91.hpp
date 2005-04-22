#ifndef _TESTODESOLVERFORLR91_HPP_
#define _TESTODESOLVERFORLR91_HPP_


#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "LR91Model.hpp"
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "SodiumCurrentLR91.hpp"
#include "PotassiumTimeDependentCurrentLR91.hpp"
#include "PotassiumTimeIndependentCurrentLR91.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"


class TestOdeSolverForLR91 : public CxxTest::TestSuite
{
    public:
    
    // Test Ode Solver for LR91
    void testOdeSolverForLR91(void)
    {
    //Tests that setting and getting stimulus works -test passed
//        double time = 0.8;
//        double magnitudeOfStimulus = 80.0;        
//        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
//        double stimulusValue = pStimulus->GetStimulus(time);
//        std::cout<<"Stimulus value is "<< stimulusValue<< std::endl;
//        
//        // test 1
       // double time = 0.2;
        //double magnitudeOfStimulus = 80.0;        
      //  AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
       // LR91OdeFun lr91OdeFun(pStimulus);
//        double voltage = -75.0;
//        double m = 0.9833;
//        double h = 0.0017;
//        double j = 0.9895;
//        double d = 0.003;
//        double f = 1;
//        double x = 0.0056;
//        double caI = 0.0002;
//          
//        std::vector<double> intialConditions;
//        intialConditions.push_back(voltage);
//        intialConditions.push_back(m);
//        intialConditions.push_back(h);
//        intialConditions.push_back(j);
//        intialConditions.push_back(d);
//        intialConditions.push_back(f);
//        intialConditions.push_back(x);
//        intialConditions.push_back(caI);
        // double t = 0.2;
        //std::vector<double> Result = lr91OdeFun.EvaluateYDerivatives(t, intialConditions);
        
        
//        DO REGRESSION TEST here....
//        double iStim = stimulusValue;
//        
//          
//        TS_ASSERT_DELTA(Result[0], (-iStim -iTotal) /cMembrane , 0.0001);
//        TS_ASSERT_DELTA(Result[1],   , 0.0001);
//        TS_ASSERT_DELTA(Result[2],  , 0.0001);
//        TS_ASSERT_DELTA(Result[3],  , 0.0001);
//        TS_ASSERT_DELTA(Result[4],  , 0.0001);
//        TS_ASSERT_DELTA(Result[5],  , 0.0001);
//        TS_ASSERT_DELTA(Result[6],  , 0.0001);
//        TS_ASSERT_DELTA(Result[7],  , 0.0001);
//          
//        TS_TRACE("IT WORKS!!!! :)) ");
        // TS_TRACE("LupRudy WORKS!!!! :)) ");
        
//        
//       //test 2 
//        voltage = -40.0;
//        m = 0.1;
//        h = 0.03;
//        j = 0.7;
//        d = 0.05;
//        f = 0.5;
//        x = 0.1;
//        caI = 0.001;
//        double magnitudeOfStimulus = 80.0;        
//        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
//       
		double voltage = -75.0;
        double m = 0.9833;
        double h = 0.0017;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002;
        double magnitudeOfStimulus = -5.1;   
          
        std::vector<double> intialConditions;
        intialConditions.push_back(voltage);
        intialConditions.push_back(m);
        intialConditions.push_back(h);
        intialConditions.push_back(j);
        intialConditions.push_back(d);
        intialConditions.push_back(f);
        intialConditions.push_back(x);
        intialConditions.push_back(caI);
		             
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus); 
        
        LR91Model *pLR91Model = new LR91Model(voltage, m, h, j, d, f, x, caI, pStimulus);
        
        // TS_TRACE("LupRudy is initiated !!!! :)) ");
        
        
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        
        double startTime = 0.0;
        double endTime = 10.0;
        double timeStep = 0.001;             
        
        OdeSolution Solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
                               intialConditions, pMySolver);
                            
                               
//       for (int i=0; i < Solution.mSolutions.size(); i++)
//       {
//       	    std::cout << Solution.mSolutions[i][0] << std::endl; 
//       	  //  sleep(1);
//       }                     
                               
        //Solution.SaveToFile("LRresult.dat");
                               
//        output to matlab file
        
        // TS_TRACE("MOREOVER OdeSolver SOLVES!!!! :)) ");

    }
    
    
    
   
	// Tests that Na gating variables within range 0 to 1 inclusive
	void testSodiumGatingVariables( void )
	{
		double m = 0.0017;
		double h = 0.9833;
		double j = 0.9895;
		double voltage = -75.0;
		
		SodiumCurrentLR91 *pINa = new SodiumCurrentLR91();
		pINa->UpdateMagnitudeOfCurrent(voltage,m,h,j);
		
		std::cout << "\n";
		std::cout << "INa m gate: " << pINa->GetM() << "\n";
		TS_ASSERT(pINa->GetM() >= 0);
		TS_ASSERT(pINa->GetM() <= 1);
		
		std::cout << "INa h gate: " << pINa->GetH() << "\n";
		TS_ASSERT(pINa->GetH() >= 0);
		TS_ASSERT(pINa->GetH() <= 1);		
	}
	
	// Tests that K gating variables within range 0 to 1 inclusive
	void testPotassiumTimeDependentGatingVariables( void )
	{
		double x = 0.0056;
		double voltage = -75.0;
		
		PotassiumTimeDependentCurrentLR91 *pIK = new PotassiumTimeDependentCurrentLR91();
		pIK->UpdateMagnitudeOfCurrent(voltage,x);
			
		std::cout << "IK n gate: " << pIK->GetX() << "\n";
		TS_ASSERT(pIK->GetX() >= 0);
		TS_ASSERT(pIK->GetX() <= 1);		
	}
	
	// Tests that K gating variables within range 0 to 1 inclusive
	void testPotassiumTimeIndependentGatingVariables( void )
	{
		double voltage = -75.0;
		
		PotassiumTimeIndependentCurrentLR91 *pIK = new PotassiumTimeIndependentCurrentLR91();
		pIK->UpdateMagnitudeOfCurrent(voltage);
			
		std::cout << "IK n gate: " << pIK->GetK1(voltage) << "\n";
		TS_ASSERT(pIK->GetK1(voltage) >= 0);
		TS_ASSERT(pIK->GetK1(voltage) <= 1);		
	}
	
	// Tests that K gating variables within range 0 to 1 inclusive
	void testPlateauPotassiumCurrentLR91( void )
	{
		double x = 0.0056;
		double voltage = -75.0;
		
		PotassiumTimeDependentCurrentLR91 *pIK = new PotassiumTimeDependentCurrentLR91();
		pIK->UpdateMagnitudeOfCurrent(voltage,x);
			
		std::cout << "IK n gate: " << pIK->GetX() << "\n";
		TS_ASSERT(pIK->GetX() >= 0);
		TS_ASSERT(pIK->GetX() <= 1);		
	}
	
	
	
    
};



#endif //_TESTODESOLVERFORLR91_HPP_
