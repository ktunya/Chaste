#ifndef CELLLUO_RUDY_1991FROMCELLMLOPT_HPP_
#define CELLLUO_RUDY_1991FROMCELLMLOPT_HPP_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: luo_rudy_1991
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translate: 8117, pycml: 8095)
//! on Fri Mar  5 16:15:37 2010
//! 
//! <autogenerated>

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"

class Cellluo_rudy_1991FromCellMLOpt : public AbstractCardiacCell
{
    // 
    // Settable parameters and readable variables
    // 
    double var_membrane__I_stim;
    double var_membrane__i_Na;
    double var_membrane__i_si;
    double var_membrane__i_K;
    double var_membrane__i_K1;
    double var_membrane__i_Kp;
    double var_membrane__i_b;
    
public:
    double Get_membrane__I_stim();
    double Get_membrane__i_Na();
    double Get_membrane__i_si();
    double Get_membrane__i_K();
    double Get_membrane__i_K1();
    double Get_membrane__i_Kp();
    double Get_membrane__i_b();
    Cellluo_rudy_1991FromCellMLOpt(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~Cellluo_rudy_1991FromCellMLOpt();
    void VerifyStateVariables();
    double GetIIonic();
    void EvaluateYDerivatives(double var_environment__time, const std::vector<double>& rY, std::vector<double>& rDY);
};


#endif // CELLLUO_RUDY_1991FROMCELLMLOPT_HPP_
