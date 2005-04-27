#ifndef _TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include "LinearHeatEquationPde.hpp"
#include "SimpleLinearSolver.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <iostream>
#include "Node.hpp" 
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"

class TestSimpleLinearEllipticAssembler : public CxxTest::TestSuite 
{
	public:
	
	void TestAssembleOnElement( void )
	{
		LinearHeatEquationPde<1> pde;
		std::vector<Node<1>*> nodes;
		nodes.push_back(new Node<1>(0, false, 1.0));
		nodes.push_back(new Node<1>(1, false, 3));
		Element<1,1> element(nodes);
		LinearBasisFunction<1> basis_function;
		MatrixDouble ael(2,2);
		VectorDouble bel(2);
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		assembler.AssembleOnElement(element, ael, bel, &pde, basis_function);
		
		TS_ASSERT_DELTA(ael(0,0),0.5, 1e-12);
		TS_ASSERT_DELTA(ael(0,1),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,1),0.5, 1e-12);
		
		TS_ASSERT_DELTA(bel(0),1, 1e-12);
		TS_ASSERT_DELTA(bel(1),1, 1e-12);
		
	}
	
	void TestAssembleOnElement2DCanonical ( void )
	{
		LinearHeatEquationPde<2> pde;
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
		nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
		nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
		Element<2,2> element(nodes);
		LinearBasisFunction<2> basis_function;
		MatrixDouble ael(3,3);
		VectorDouble bel(3);
		
		SimpleLinearEllipticAssembler<2,2> assembler;
		assembler.AssembleOnElement(element, ael, bel, &pde, basis_function);
		
		TS_ASSERT_DELTA(ael(0,0),1.0, 1e-12);
		TS_ASSERT_DELTA(ael(0,1),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(0,2),-0.5, 1e-12);
		
		TS_ASSERT_DELTA(ael(1,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,1),0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,2),0.0, 1e-12);
		
		
		TS_ASSERT_DELTA(ael(2,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(2,1),0.0, 1e-12);
		TS_ASSERT_DELTA(ael(2,2),0.5, 1e-12);
		
		TS_ASSERT_DELTA(bel(0),1.0/6.0, 1e-12);
		TS_ASSERT_DELTA(bel(1),1.0/6.0, 1e-12);
		TS_ASSERT_DELTA(bel(2),1.0/6.0, 1e-12);
		
	}
    
	void TestAssembleOnElement2DGeneral ( void )
	{
		LinearHeatEquationPde<2> pde;
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
		nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
		nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
		Element<2,2> element(nodes);
		LinearBasisFunction<2> basis_function;
		MatrixDouble ael(3,3);
		VectorDouble bel(3);
		
		SimpleLinearEllipticAssembler<2,2> assembler;
		assembler.AssembleOnElement(element, ael, bel, &pde, basis_function);
		
		TS_ASSERT_DELTA(ael(0,0),1.0, 1e-12);
		TS_ASSERT_DELTA(ael(0,1),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(0,2),-0.5, 1e-12);
		
		TS_ASSERT_DELTA(ael(1,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,1),0.5, 1e-12);
		TS_ASSERT_DELTA(ael(1,2),0.0, 1e-12);
		
		
		TS_ASSERT_DELTA(ael(2,0),-0.5, 1e-12);
		TS_ASSERT_DELTA(ael(2,1),0.0, 1e-12);
		TS_ASSERT_DELTA(ael(2,2),0.5, 1e-12);
		
		TS_ASSERT_DELTA(bel(0),5.0/6.0, 1e-12);
		TS_ASSERT_DELTA(bel(1),5.0/6.0, 1e-12);
		TS_ASSERT_DELTA(bel(2),5.0/6.0, 1e-12);
		
	}
	
	void TestWithHeatEquation()   
	{ 
		PetscInitialize(0, NULL, 0, 0);
		
		// Create mesh (by hand!)
		const int num_elements = 10;
		ConformingTetrahedralMesh<1,1> mesh(num_elements);
		std::vector<Node<1>*> nodes;
		for (int i=0; i<num_elements+1; i++)
		{
			nodes.push_back(new Node<1>(i, true, 0.0 + 0.15*i));
			mesh.AddNode(*nodes[i]);
		}
		for (int i=0; i<num_elements; i++)
		{
			std::vector<Node<1>*> element_nodes;
			element_nodes.push_back(nodes[i]);
			element_nodes.push_back(nodes[i+1]);
			Element<1,1> element(element_nodes);
			mesh.AddElement(element);
		}
		
		// Instantiate PDE object
		LinearHeatEquationPde<1> pde;  
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution should be u = 0.5*x*(3-x)
		for (int i=0; i < num_elements+1; i++)
		{
			double x = 0.0 + 0.15*i;
			double u = 0.5*x*(3-x);
			TS_ASSERT_DELTA(res[i], u, 0.001);
		}
		VecRestoreArray(result, &res);
	}
    
	void TestWithHeatEquationAndMeshReader()   
	{ 
		PetscInitialize(0, NULL, 0, 0);
		
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/trivial_1d_mesh");
		ConformingTetrahedralMesh<1,1> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		LinearHeatEquationPde<1> pde;  
		
		// Boundary conditions
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<1,1> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution should be u = 0.5*x*(3-x)
		for (int i=0; i < mesh.GetNumElements()+1; i++)
		{
			double x = 0.0 + 0.15*i;
			double u = 0.5*x*(3-x);
			TS_ASSERT_DELTA(res[i], u, 0.001);
		}
		VecRestoreArray(result, &res);
	}

    void TestWithHeatEquation2()
    {
        PetscInitialize(0, NULL, 0, 0);
        
        // Create mesh (by hand!)
        const int num_elements = 5;
        ConformingTetrahedralMesh<1,1> mesh(num_elements);
        
        Node<1> *right_hand_node = new Node<1>(0, true, -1.0);
        mesh.AddNode(*right_hand_node);
        
        
        for (int i=1; i<num_elements; i++)
        {
            std::vector<Node<1>*> element_nodes;
            Node<1> *left_hand_node = new Node<1>(i, false, -1.0-(2.0*i/num_elements));
            mesh.AddNode(*left_hand_node);
            element_nodes.push_back(left_hand_node);
            element_nodes.push_back(right_hand_node);
            Element<1,1> element(element_nodes);
            mesh.AddElement(element);
            right_hand_node = left_hand_node;
        }
        
		std::vector<Node<1>*> element_nodes;
		Node<1> *left_hand_node = new Node<1>(num_elements, true, -3.0);
        mesh.AddNode(*left_hand_node);
        element_nodes.push_back(left_hand_node);
        element_nodes.push_back(right_hand_node);
        Element<1,1> element(element_nodes, true);
        mesh.AddElement(element);
        
        // Instantiate PDE object
        LinearHeatEquationPde<1> pde;
        
        // Boundary conditions
        // u(-1)=1 u'(-3)=0
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        //std::cout << *mesh.GetNodeAt(0)<< "\n";
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        
        ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(0.0);
        // element is still the last element in the mesh
        mesh.AddSurfaceElement(element.GetLowerOrderElement(1)); 
        bcc.AddNeumannBoundaryCondition(element.GetLowerOrderElement(1),
                                        pNeumannBoundaryCondition);
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler 
        SimpleLinearEllipticAssembler<1,1> assembler;   
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
        
        
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int i=0; i<=num_elements; i++)//result.Size()
        {
            double x = mesh.GetNodeAt(i)->GetPoint()[0];
        }
        for (int i=0; i < num_elements+1; i++)
        {
            double x = -1.0- 0.4*i;
            double u = 1 - 0.5*(x+1)*(5+x);
            TS_ASSERT_DELTA(res[i], u, 0.001);
        }
    }
    
    
    
    void TestWithHeatEquationNonzeroNeumannCondition()
    {
        PetscInitialize(0, NULL, 0, 0);
        
        // Create mesh (by hand!)
        const int num_elements = 5;
        ConformingTetrahedralMesh<1,1> mesh(num_elements);
        
        Node<1> *right_hand_node = new Node<1>(0, true, -1.0);
        mesh.AddNode(*right_hand_node);
        
        
        for (int i=1; i<num_elements; i++)
        {
            std::vector<Node<1>*> element_nodes;
            Node<1> *left_hand_node = new Node<1>(i, false, -1.0-(2.0*i/num_elements));
            mesh.AddNode(*left_hand_node);
            element_nodes.push_back(left_hand_node);
            element_nodes.push_back(right_hand_node);
            Element<1,1> element(element_nodes);
            mesh.AddElement(element);
            right_hand_node = left_hand_node;
        }
        std::vector<Node<1>*> element_nodes;
        Node<1> *left_hand_node = new Node<1>(num_elements, true, -3.0);
        mesh.AddNode(*left_hand_node);
        element_nodes.push_back(left_hand_node);
        element_nodes.push_back(right_hand_node);
        Element<1,1> element(element_nodes, true);
        mesh.AddElement(element);
        
        // Instantiate PDE object
        LinearHeatEquationPde<1> pde;
        
        // Boundary conditions
        // u(-1)=1 u'(-3)=1
        BoundaryConditionsContainer<1,1> bcc;
        ConstBoundaryCondition<1>* pBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[0], -1, 1e-12);
        
        ConstBoundaryCondition<1>* pNeumannBoundaryCondition = new ConstBoundaryCondition<1>(1.0);
        // element is still the last element in the mesh
        mesh.AddSurfaceElement(element.GetLowerOrderElement(1));
        bcc.AddNeumannBoundaryCondition(element.GetLowerOrderElement(1),
                                        pNeumannBoundaryCondition);
        TS_ASSERT_DELTA(element.GetLowerOrderElement(1)->GetNodeLocation(0,0), -3, 1e-12);
        
        // Linear solver
        SimpleLinearSolver solver;
        
        // Assembler 
        SimpleLinearEllipticAssembler<1,1> assembler;
        
        Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);       
        
        double *res;
        int ierr = VecGetArray(result, &res);
        for (int i=0; i<=num_elements; i++)
        {
            double x = mesh.GetNodeAt(i)->GetPoint()[0];
        }
        for (int i=0; i < num_elements+1; i++)
        {
            double x = -1.0- 0.4*i;
            double u = -0.5*x*x - 2*x - 0.5;
            TS_ASSERT_DELTA(res[i], u, 0.001);
        }
        //TS_TRACE("here simp lin");
    }

	/**
	 * \todo
	 * Skeleton 2d test. Need to check real solution.
	 */
	void Test2dHeatEquation()
	{
		PetscInitialize(0, NULL, 0, 0);
		
		// Create mesh from mesh reader
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/square_4_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Instantiate PDE object
		LinearHeatEquationPde<2> pde;
		
		// Boundary conditions
        BoundaryConditionsContainer<2,2> bcc;
        ConstBoundaryCondition<2>* pBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(0), pBoundaryCondition);
        pBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), pBoundaryCondition);
        pBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(2), pBoundaryCondition);
        pBoundaryCondition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNodeAt(3), pBoundaryCondition);
        
		// Linear solver
		SimpleLinearSolver solver;
		
		// Assembler
		SimpleLinearEllipticAssembler<2,2> assembler;
		
		Vec result = assembler.AssembleSystem(mesh, &pde, bcc, &solver);
		
		// Check result
		double *res;
		int ierr = VecGetArray(result, &res);
		// Solution at 4th node should be 0.08
		TS_ASSERT_DELTA(res[4], 1.0/12.0, 0.001);
	}
	

};
 
#endif //_TESTSIMPLELINEARELLIPTICASSEMBLER_HPP_
