#ifndef TESTPERFORMANCE_HPP_
#define TESTPERFORMANCE_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "PerformanceTester.hpp"

class TestPerformance : public CxxTest::TestSuite
{   
public:

    
    void TestPerf() throw(Exception)
    {
        // write headings
        PerformanceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2>::DisplayHeadings();
        EventHandler::Headings();

        // base line
        PerformanceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2> tester;
        tester.SimTime=0.0;
        tester.MeshNum=3;
        tester.Run();
        
        EventHandler::Report();
    }
};

#endif /*TESTPERFORMANCE_HPP_*/
