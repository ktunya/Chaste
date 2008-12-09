/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTVERTEXMESH_HPP_
#define TESTVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "VertexElement.hpp"
#include "VertexMesh.hpp"

class TestVertexMesh : public CxxTest::TestSuite
{
public:
    void TestBasicVertexMesh()
    {
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        basic_nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        
        std::vector<Node<2>*> nodes_elem_1, nodes_elem_2;
        
        nodes_elem_1.push_back(basic_nodes[0]);
        nodes_elem_1.push_back(basic_nodes[1]);
        nodes_elem_1.push_back(basic_nodes[2]);
        
        nodes_elem_2.push_back(basic_nodes[0]);
        nodes_elem_2.push_back(basic_nodes[2]);
        nodes_elem_2.push_back(basic_nodes[3]);
        
        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        basic_vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));
        
        VertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);
        
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumVertexElements(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 4u);
                       
        std::set<unsigned> temp_list1;
        
        temp_list1.insert(0u);
        
        TS_ASSERT_EQUALS(basic_nodes[1]->rGetContainingElementIndices(), temp_list1);
        
        temp_list1.insert(1u);
        
        TS_ASSERT_EQUALS(basic_nodes[0]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(basic_nodes[2]->rGetContainingElementIndices(), temp_list1);
        
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(basic_nodes[3]->rGetContainingElementIndices(), temp_list2);
        
//        temp_vertex_elements_owned_by_nodes1.push_back(basic_vertex_elements[0]);
//        temp_vertex_elements_owned_by_nodes1.push_back(basic_vertex_elements[1]);
//        
//        temp_vertex_elements_owned_by_nodes2.push_back(basic_vertex_elements[0]);
//        
//        temp_vertex_elements_owned_by_nodes3.push_back(basic_vertex_elements[1]);
//               
//        vertex_elements_owned_by_nodes.push_back(temp_vertex_elements_owned_by_nodes1);
//        vertex_elements_owned_by_nodes.push_back(temp_vertex_elements_owned_by_nodes2);
//        vertex_elements_owned_by_nodes.push_back(temp_vertex_elements_owned_by_nodes1);
//        vertex_elements_owned_by_nodes.push_back(temp_vertex_elements_owned_by_nodes3);
//        
//        std::set<unsigned> vertex_elements_owned_by_node_returned;
//        
//        vertex_elements_owned_by_node_returned = basic_nodes[0]->rGetContainingElementIndices;
//        
//        TS_ASSERT_EQUALS(vertex_elements_owned_by_node_returned[0], basic_vertex_elements[0]);
//        std::cout<< "elem number " << vertex_elements_owned_by_node_returned[1]->GetIndex() << std::endl;
//        TS_ASSERT_EQUALS(vertex_elements_owned_by_node_returned[1], basic_vertex_elements[1]);
//        
//        vertex_elements_owned_by_node_returned = basic_vertex_mesh.GetElementsOwnedByNode(basic_nodes[1]);
//        
//        TS_ASSERT_EQUALS(vertex_elements_owned_by_node_returned[0], basic_vertex_elements[0]);
//        
//        vertex_elements_owned_by_node_returned = basic_vertex_mesh.GetElementsOwnedByNode(basic_nodes[2]);
//        
//        TS_ASSERT_EQUALS(vertex_elements_owned_by_node_returned[0], basic_vertex_elements[0]);
//        TS_ASSERT_EQUALS(vertex_elements_owned_by_node_returned[1], basic_vertex_elements[1]);
//        
//        vertex_elements_owned_by_node_returned = basic_vertex_mesh.GetElementsOwnedByNode(basic_nodes[3]);
//        
//        TS_ASSERT_EQUALS(vertex_elements_owned_by_node_returned[0], basic_vertex_elements[1]);
    }
};    

#endif /*TESTVERTEXMESH_HPP_*/
