/**
 * Concrete TestOde1 class
 */ 
#ifndef _TESTODE1_HPP_
#define _TESTODE1_HPP_
// TestOde1.hpp
#include "AbstractOdeSystem.hpp"


class TestOde1 : public AbstractOdeSystem
{
	public :
	TestOde1();
	
	std::vector<double> TestOde1::EvaluateYDerivatives (double time, const std::vector<double> &rY);

};

#endif //_TESTODE1_HPP_
