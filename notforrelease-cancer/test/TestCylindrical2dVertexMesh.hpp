/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTCYLINDRICAL2DVERTEXMESH_HPP_
#define TESTCYLINDRICAL2DVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Cylindrical2dVertexMesh.hpp"


class TestCylindrical2dVertexMesh : public CxxTest::TestSuite
{
public:

    void TestMeshGetWidth(void)
    {
        // Create mesh
        Cylindrical2dVertexMesh cylindrical_vertex_mesh(3, 3, 0.01, 2.0);

        // Test GetWidthExtremes() method
        c_vector<double,2> width_extremes = cylindrical_vertex_mesh.GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = cylindrical_vertex_mesh.GetWidthExtremes(1u);

        TS_ASSERT_DELTA(width_extremes[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(width_extremes[1], 2.8867, 1e-4);

        TS_ASSERT_DELTA(height_extremes[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(height_extremes[1], 3.5000, 1e-4);

        // Test GetWidth() method
        double width = cylindrical_vertex_mesh.GetWidth(0);
        double height = cylindrical_vertex_mesh.GetWidth(1);

        TS_ASSERT_DELTA(width, 2.8867, 1e-4);
        TS_ASSERT_DELTA(height, 3.5000, 1e-4);
    }

    void TestGetVectorFromAtoB() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(6, 6, 0.01, 2.0);

        c_vector<double, 2> node36_location = mesh.GetNode(36)->rGetLocation();
        c_vector<double, 2> node38_location = mesh.GetNode(38)->rGetLocation();

        // Test a normal vector and distance calculation
        c_vector<double, 2> vector = mesh.GetVectorFromAtoB(node36_location, node38_location);
        TS_ASSERT_DELTA(vector[0], 1.7320, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1.7320, 1e-4);
        TS_ASSERT_DELTA(mesh.GetDistanceBetweenNodes(36, 38), 1.7320, 1e-4);

        // Test the opposite vector
        vector = mesh.GetVectorFromAtoB(node38_location, node36_location);
        TS_ASSERT_DELTA(vector[0], -1.7320, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);

        // Test a periodic calculation

        c_vector<double, 2> node25_location = mesh.GetNode(25)->rGetLocation();
        c_vector<double, 2> node28_location = mesh.GetNode(28)->rGetLocation();

        vector = mesh.GetVectorFromAtoB(node25_location, node28_location);
        TS_ASSERT_DELTA(vector[0], 1.7320, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.5000, 1e-4);

        /// \todo add more cases - see also TestCylindrical2dMesh.hpp (#918)
    }


    void TestSetNodeLocationForCylindricalMesh() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(6, 6, 0.01, 2.0);

        // Move one of the nodes to near the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = -0.01;
        new_point_location[1] = 2.5;
        ChastePoint<2> new_point(new_point_location);

        // This node was on left and is now near the right
        mesh.SetNode(34u, new_point);
        TS_ASSERT_DELTA(mesh.GetNode(34u)->rGetLocation()[0], 5.4748, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(34u)->rGetLocation()[1], 2.5, 1e-4);

        /// \todo add more cases - see also TestCylindrical2dMesh.hpp (#918)
    }


    void TestAddNodeAndReMesh() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(6, 6, 0.01, 2.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 96u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 36u);

        // Choose a node on the left boundary
        ChastePoint<2> point = mesh.GetNode(62)->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(point[1], 4.5000, 1e-4);

        // Create a new node close to this node
        point.SetCoordinate(0, -0.01);
        point.SetCoordinate(1, 4.5);
        Node<2>* p_node = new Node<2>(mesh.GetNumNodes(), point);

        unsigned old_num_nodes = mesh.GetNumNodes();

        // Add this new node to the mesh
        unsigned new_index = mesh.AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Remesh to update correspondences
        mesh.SetEdgeDivisionThreshold(1000); // set high threshold to avoid more nodes appearing in the remesh
        NodeMap map(mesh.GetNumElements());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 97u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 36u);

        TS_ASSERT_DELTA(mesh.GetNode(new_index)->rGetLocation()[0], 5.4748, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(new_index)->rGetLocation()[1], 4.5000, 1e-4);

        // Now tet AddNode() when mDeletedNodeIndices is populated

        // Label node 29 as deleted
        mesh.mDeletedNodeIndices.push_back(29);

        // Create a new node close to this node
        ChastePoint<2> point2;
        point2.SetCoordinate(0, 2.0);
        point2.SetCoordinate(1, 2.1);
        Node<2>* p_node2 = new Node<2>(mesh.GetNumNodes(), point);

        // Add this new node to the mesh
        new_index = mesh.AddNode(p_node2);
        TS_ASSERT_EQUALS(new_index, 29u);

        /// \todo add more cases - see also TestCylindrical2dMesh.hpp (#918)
    }

    /// \todo add archiving test - see also TestCylindrical2dMesh.hpp (#918)

};

#endif /*TESTCYLINDRICAL2DVERTEXMESH_HPP_*/
