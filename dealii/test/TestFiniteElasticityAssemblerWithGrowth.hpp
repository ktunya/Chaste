#ifndef TESTFINITEELASTICITYASSEMBLERWITHGROWTH_HPP_
#define TESTFINITEELASTICITYASSEMBLERWITHGROWTH_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include "FiniteElasticityAssemblerWithGrowth.cpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

#include "MooneyRivlinMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "ExponentialMaterialLaw.hpp"

#include "FiniteElasticityTools.hpp"
#include "ConcentrationBasedTumourSourceModel.hpp"
#include "ConstantTumourSourceModel.hpp"

#include "grid/tria_boundary_lib.h"



// todos: proper test of answers, compare numerical jacobian
// sensible test once s set up


class TestFiniteElasticityAssemblerWithGrowth : public CxxTest::TestSuite
{
public :
    void TestExceptions() throw(Exception)
    {
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        ConstantTumourSourceModel<2> source_model(1.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        // elements haven't been set as GROWING_REGION or NON_GROWING_REGION
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssemblerWithGrowth<2> bad_fe_with_growth(&mesh,&mooney_rivlin_law,body_force,1.0,"",&source_model));
        
        FiniteElasticityTools<2>::SetAllElementsAsNonGrowingRegion(mesh);
        
        // no elements set as GROWING_REGION
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssemblerWithGrowth<2> bad_fe_with_growth2(&mesh,&mooney_rivlin_law,body_force,1.0,"",&source_model));
        
        // set the first element as growing
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        element_iter->set_material_id(GROWING_REGION);
        
        // should construct ok now
        FiniteElasticityAssemblerWithGrowth<2> fe_with_growth(&mesh,&mooney_rivlin_law,body_force,1.0,"",&source_model);
        
        // set times not been called
        TS_ASSERT_THROWS_ANYTHING(fe_with_growth.Run());
        
        // start time > end time
        TS_ASSERT_THROWS_ANYTHING(fe_with_growth.SetTimes(1.0, 0.0, 0.01));
        
        // dt negative
        TS_ASSERT_THROWS_ANYTHING(fe_with_growth.SetTimes(0.0, 1.0, -0.01));
        
        // none of the above should throw now
        TS_ASSERT_THROWS_NOTHING(fe_with_growth.SetTimes(0.0, 1.0, 0.01));
    }
    
    void TestWithSimpleProblem() throw(Exception)
    {
        Vector<double> body_force(2); // zero
        double density = 1.233;
        
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);
        
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        double initial_elem_volume = 1.0/mesh.n_active_cells();
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        // set all elements as growing
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, zero, 100);
        double source_value = 2;        
        ConstantTumourSourceModel<2> source_model(source_value);
        
        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/simple",
                                                                      &source_model);
                
                
        // loop over all the elements, and if it is in the growing region, check
        // each node has an ode system associated with it...
        TriangulationVertexIterator<2> vertex_iter(&mesh);
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            TS_ASSERT_EQUALS(finiteelas_with_growth.IsGrowingNode(vertex_index), true);
            vertex_iter.Next();
        }
        
        // run
        double end_time = 0.2;
        finiteelas_with_growth.SetTimes(0.0, end_time, 0.1);
        finiteelas_with_growth.Run();
 
        std::vector<Vector<double> >& deformed_position
            = finiteelas_with_growth.rGetDeformedPosition();
 
        // test
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            double x_undef_node_0 = element_iter->vertex(0)[0];
            double y_undef_node_0 = element_iter->vertex(0)[1];
            
            // look at elements in the top right-hand quadrant, ie those away
            // from the fixed corner (0,0). (if there was no fixed nodes, then
            // the solution would just be iostropic enlargement (simple stretching)  
            if( (x_undef_node_0>0.5) && (y_undef_node_0>0.5))
            {
                double x0 = deformed_position[0](element_iter->vertex_index(0));
                double y0 = deformed_position[1](element_iter->vertex_index(0));

                double x1 = deformed_position[0](element_iter->vertex_index(1));
                double y1 = deformed_position[1](element_iter->vertex_index(1));

                double x2 = deformed_position[0](element_iter->vertex_index(2));
                double y2 = deformed_position[1](element_iter->vertex_index(2));

                double x3 = deformed_position[0](element_iter->vertex_index(3));
                double y3 = deformed_position[1](element_iter->vertex_index(3));

                // these elements are away from the fixed corner, so should be very 
                // like enlarged rectangles. Verify this, ie check x0=x3, x1=x2, etc
                TS_ASSERT_DELTA(x0, x3, 1e-3);
                TS_ASSERT_DELTA(x1, x2, 1e-3);
                TS_ASSERT_DELTA(y0, y1, 1e-3);
                TS_ASSERT_DELTA(y2, y3, 1e-3);

                // check volume of enlarged square is as expected
                // dgdt = 0.5 g rho s, so g = exp(0.5 rho s T)
                // F = gI so detF = g^2 = exp(rho s T)
                double expected_volume = exp(density*end_time*source_value)*initial_elem_volume;
                TS_ASSERT_DELTA( (x1-x0)*(y2-y1), expected_volume, 1e-3);
            }
            element_iter++;
        }
    }
};
#endif /*TESTFINITEELASTICITYASSEMBLERWITHGROWTH_HPP_*/
