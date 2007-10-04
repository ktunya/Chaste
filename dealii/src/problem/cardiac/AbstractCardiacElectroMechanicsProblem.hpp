#ifndef ABSTRACTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define ABSTRACTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include "MonodomainProblem.hpp"
#include "AbstractCardiacMechanicsAssembler.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "NhsCellularMechanicsOdeSystem.hpp"
#include "FiniteElasticityTools.hpp"
#include "AbstractElasticityAssembler.hpp"

/* todos:
 * 
 * add comments
 * add tests
 * 
 * det F
 * 
 * move mesh stuff out
 * 
 * think about architecture (of AbstractCardiacProblem) when this is done properly..
 */



/**
 *  At the beginning of a two mesh simulation we need to figure out and store
 *  which (electrics-mesh) element each (mechanics-mesh) gauss point is in, and
 *  what the weight of that gauss point for that particular element is. This struct
 *  just contains this two pieces of data
 */
template<unsigned DIM>
struct ElementAndWeights
{
    unsigned ElementNum;
    c_vector<double,DIM+1> Weights;  
};


template<unsigned DIM>
class AbstractCardiacElectroMechanicsProblem
{
protected :
    /*< The cardiac problem class */
    MonodomainProblem<DIM>* mpMonodomainProblem;
    
    /*< The mechanics assembler */
    AbstractCardiacMechanicsAssembler<DIM>* mpCardiacMechAssembler;  

    /*< End time. The start time is assumed to be 0.0 */
    double mEndTime;
    /*< The timestep. TODO: different timesteps for different bits */
    double mTimeStep;    
    
    /*< A chaste mesh for the electrics */
    ConformingTetrahedralMesh<DIM,DIM>* mpElectricsMesh;
    /*<  A dealii mesh for the mechanics */
    Triangulation<DIM>*                 mpMechanicsMesh;

    /** 
     *  The (electrics-mesh) element numbers and saying which element each 
     *  (mechanics-mesh) gauss point is in, and the weight of that gauss point 
     *  for that particular element is.
     */
    std::vector<ElementAndWeights<DIM> > mElementAndWeightsForQuadPoints;

    /*< Whether to use an explicit or implicit method */
    bool mUseExplicitMethod;

    /*< Output directory, relative to TEST_OUTPUT */
    std::string mOutputDirectory;
    /*< Whether to write any output */
    bool mWriteOutput;
    /*< when to write output */    
    const static int WRITE_EVERY_NTH_TIME = 1; 
    
    virtual void ConstructMechanicsAssembler()=0;
    virtual void ConstructMeshes()=0;

public :
    AbstractCardiacElectroMechanicsProblem(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                           double endTime,
                                           double timeStep,
                                           bool useExplicitMethod,
                                           std::string outputDirectory = "")
    {
        assert(pCellFactory != NULL);
        mpMonodomainProblem = new MonodomainProblem<DIM>(pCellFactory);
        
        assert(endTime > 0);
        mEndTime = endTime;
        mTimeStep = timeStep;
        
        // check whether output is required
        mWriteOutput = (outputDirectory!="");
        mOutputDirectory = outputDirectory;
                
        mUseExplicitMethod = useExplicitMethod;
        mpCardiacMechAssembler = NULL;
    }   
    
    void Initialise()
    {
        ConstructMeshes();     
                        
        mpMonodomainProblem->SetMesh(mpElectricsMesh);
        mpMonodomainProblem->Initialise();

        ConstructMechanicsAssembler();

        // find the element nums and weights for each gauss point in the mechanics mesh
        mElementAndWeightsForQuadPoints.resize(mpCardiacMechAssembler->GetTotalNumQuadPoints());

        std::vector<std::vector<double> > quad_point_posns
           = FiniteElasticityTools<DIM>::GetQuadPointPositions(*mpMechanicsMesh, mpCardiacMechAssembler->GetNumQuadPointsInEachDimension());

        
        for(unsigned i=0; i<quad_point_posns.size(); i++)
        {
            ChastePoint<DIM> point;

            for(unsigned j=0;j<DIM;j++)
            {
                point.rGetLocation()[j]=quad_point_posns[i][j];
                std::cout << point[j] << " ";
            }
            
            std::cout << "\n";
            
            unsigned elem_index = mpElectricsMesh->GetContainingElementIndex(point);
            c_vector<double,DIM+1> weight = mpElectricsMesh->GetElement(elem_index)->CalculateInterpolationWeights(point);
            
            mElementAndWeightsForQuadPoints[i].ElementNum = elem_index;
            mElementAndWeightsForQuadPoints[i].Weights = weight;
        }
    }

    /** 
     *  Solve the electromechanincs problem
     */    
    void Solve()
    {
        if(mpCardiacMechAssembler==NULL)
        {
            Initialise();
        }
        
        // get an electrics assembler from the problem. Note that we don't call
        // Solve() on the CardiacProblem class, we do the looping here.
        AbstractDynamicAssemblerMixin<DIM,DIM,1>* mpElectricsAssembler 
           = mpMonodomainProblem->CreateAssembler();

        // set up initial voltage etc
        Vec voltage;        
        Vec initial_voltage = mpMonodomainProblem->CreateInitialCondition();

        // create stores of lambda, lambda_dot and old lambda
        unsigned num_quad_points = mpCardiacMechAssembler->GetTotalNumQuadPoints();

        // these are only needed if explicit
        std::vector<double> lambda;
        std::vector<double> old_lambda;
        std::vector<double> dlambda_dt;
        std::vector<NhsCellularMechanicsOdeSystem> cellmech_systems;
        EulerIvpOdeSolver euler_solver;

        // this is the active tension if explicit and the calcium conc if implicit
        std::vector<double> forcing_quantity(num_quad_points,0.0);
        
        // initial cellmechanics systems, lambda, etc, if required
        if(mUseExplicitMethod)
        {
            lambda.resize(num_quad_points, 1.0);
            old_lambda.resize(num_quad_points, 1.0);
            dlambda_dt.resize(num_quad_points, 0.0);
            cellmech_systems.resize(num_quad_points);
        }

        unsigned mech_writer_counter = 0;

        // write initial positions
        if(mWriteOutput)
        {
            OutputFileHandler output_file_handler(mOutputDirectory, true);
            out_stream p_file = output_file_handler.OpenOutputFile("results_", mech_writer_counter, ".dat");
            std::vector<Vector<double> >& deformed_position = dynamic_cast<AbstractElasticityAssembler<DIM>*>(mpCardiacMechAssembler)->rGetDeformedPosition();
            for(unsigned i=0; i<deformed_position[0].size(); i++)
            {
                for(unsigned j=0; j<DIM; j++)
                {
                    (*p_file) << deformed_position[j](i) << " ";
                }
                (*p_file) << "\n";
            }
        }


        unsigned counter = 0;

        TimeStepper stepper(0.0, mEndTime, mTimeStep);
        while ( !stepper.IsTimeAtEnd() )
        {
            std::cout << "**Time = " << stepper.GetTime() << "\n" << std::flush;
            
            // solve the electrics
            mpElectricsAssembler->SetTimes(stepper.GetTime(), stepper.GetNextTime(), mTimeStep);
            mpElectricsAssembler->SetInitialCondition( initial_voltage );
            voltage = mpElectricsAssembler->Solve();

            VecDestroy(initial_voltage);
            initial_voltage = voltage;
            
            // compute Ca_I at each quad point (by interpolation, using the info on which
            // electrics element the quad point is in. Then: 
            //   Explicit: Set Ca_I on the nhs systems and solve them to get the active tension
            //   Implicit: Set Ca_I on the mechanics solver
            for(unsigned i=0; i<mElementAndWeightsForQuadPoints.size(); i++)
            {
                double interpolated_Ca_I = 0;

                Element<DIM,DIM>& element = *(mpElectricsMesh->GetElement(mElementAndWeightsForQuadPoints[i].ElementNum));
                for(unsigned node_index = 0; node_index<element.GetNumNodes(); node_index++)
                {
                    unsigned global_node_index = element.GetNodeGlobalIndex(node_index);
                    double Ca_I_at_node = mpMonodomainProblem->GetPde()->GetCardiacCell(global_node_index)->GetIntracellularCalciumConcentration();
                    interpolated_Ca_I += Ca_I_at_node*mElementAndWeightsForQuadPoints[i].Weights(node_index);
                }

                if(mUseExplicitMethod)
                {
                    // explicit: forcing quantity on the assembler is the active tension
                    cellmech_systems[i].SetLambdaAndDerivative(lambda[i], dlambda_dt[i]);
                    cellmech_systems[i].SetIntracellularCalciumConcentration(interpolated_Ca_I);
                    euler_solver.SolveAndUpdateStateVariable(&cellmech_systems[i], stepper.GetTime(), stepper.GetNextTime(), mTimeStep);
                    forcing_quantity[i] = cellmech_systems[i].GetActiveTension();
                }
                else
                {
                    // explicit: forcing quantity on the assembler is the calcium concentration
                    forcing_quantity[i] = interpolated_Ca_I;
                }
            }

            // NOTE: HERE WE SHOULD REALLY CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
            // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS.. (esp for implicit)
            
            // set the active tensions
            mpCardiacMechAssembler->SetForcingQuantity(forcing_quantity);

            // solve the mechanics
            mpCardiacMechAssembler->Solve(stepper.GetTime(), stepper.GetNextTime(), stepper.GetNextTime()-stepper.GetTime());

            // if explicit store the new lambda and update lam
            if(mUseExplicitMethod)
            {
                // update lambda and dlambda_dt;
                old_lambda = lambda;
                lambda = mpCardiacMechAssembler->rGetLambda();
                for(unsigned i=0; i<dlambda_dt.size(); i++)
                {
                    dlambda_dt[i] = (lambda[i] - old_lambda[i])/mTimeStep;
                }
            }

            // write
            if(mWriteOutput && (counter++)%WRITE_EVERY_NTH_TIME==0)
            {            
                OutputFileHandler output_file_handler(mOutputDirectory, false);
                out_stream p_file = output_file_handler.OpenOutputFile("results_", mech_writer_counter, ".dat");
                std::vector<Vector<double> >& deformed_position = dynamic_cast<AbstractElasticityAssembler<DIM>*>(mpCardiacMechAssembler)->rGetDeformedPosition();
                for(unsigned i=0; i<deformed_position[0].size(); i++)
                {
                    for(unsigned j=0; j<DIM; j++)
                    {
                        (*p_file) << deformed_position[0](i) << " ";
                    }
                    (*p_file) << "\n";
                }
                mech_writer_counter++;
            }                        
            
            // update the current time
            stepper.AdvanceOneTimeStep();
        }
        
        delete mpElectricsAssembler;
    }
};



#endif /*ABSTRACTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
