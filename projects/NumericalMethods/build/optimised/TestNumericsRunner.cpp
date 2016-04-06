/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

#include "CommandLineArguments.hpp"
int main( int argc, char *argv[] ) {
 CommandLineArguments::Instance()->p_argc = &argc;
 CommandLineArguments::Instance()->p_argv = &argv;
 return CxxTest::ErrorPrinter().run();
}
#include "projects/NumericalMethods/test/TestNumerics.hpp"

static TestNumerics suite_TestNumerics;

static CxxTest::List Tests_TestNumerics = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestNumerics( "projects/NumericalMethods/test/TestNumerics.hpp", 82, "TestNumerics", suite_TestNumerics, Tests_TestNumerics );

static class TestDescription_TestNumerics_TestWithPositionRecording : public CxxTest::RealTestDescription {
public:
 TestDescription_TestNumerics_TestWithPositionRecording() : CxxTest::RealTestDescription( Tests_TestNumerics, suiteDescription_TestNumerics, 224, "TestWithPositionRecording" ) {}
 void runTest() { suite_TestNumerics.TestWithPositionRecording(); }
} testDescription_TestNumerics_TestWithPositionRecording;

#include <cxxtest/Root.cpp>
