#ifndef TESTSPIRALPARAMETERSREADER_HPP_
#define TESTSPIRALPARAMETERSREADER_HPP_

#include <cxxtest/TestSuite.h>
#include <memory>
#include "SpiralParameters.hpp"


using std::auto_ptr;

class TestSpiralParametersReader : public CxxTest::TestSuite
{
public:
    void TestRead()
    {
        try
        {
            auto_ptr<SpiralParameters::type> p (SpiralParameters("io/test/data/Baseline.xml"));
            TS_ASSERT_EQUALS(p->SimulationDuration(), 10.0);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }
};
#endif /*TESTSPIRALPARAMETERSREADER_HPP_*/
