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
#ifndef VERTEXMESH_HPP_
#define VERTEXMESH_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>

#include "AbstractMesh.hpp"
#include "VertexElement.hpp"
#include "BoundaryElement.hpp"
#include "NodeMap.hpp"
#include "Node.hpp"
#include "Exception.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MutableMesh.hpp"
#include "VoronoiTessellation.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based tissue simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh
{
private:

    /** Vector of pointers to Nodes. */
    std::vector<Node<SPACE_DIM>*> mNodes;

    /** Vector of pointers to VertexElements. */    
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;
    
    /* Whether or not memory has been allocated for the mesh (used by destructor). */
    bool mAllocatedMemory;
    
    /** Create correspondences between VertexElements and Nodes in the mesh. */
    void SetupVertexElementsOwnedByNodes();
    
public:

    /**
     * Default constructor.
     * 
     * @param nodes vector of pointers to nodes
     * @param vertexElements  vector of pointers to VertexElements
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements);
    
    /**
     * Destructor.
     */
    ~VertexMesh();

    /**
     * Helper constructor, creates a rectangular vertex-based mesh.
     * 
     * @param numAcross number of VertexElements across
     * @param numUp number of VertexElements up
     */
    VertexMesh(unsigned numAcross, unsigned numUp);

    /**
     * @return the number of Nodes in the mesh.
     */
    unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    unsigned GetNumElements() const;

    /**
     * @param index the global index of a specified node
     * 
     * @return a pointer to the node
     */
    Node<SPACE_DIM>* GetNode(unsigned index) const;
    
    /**
     * @param index  the global index of a specified vertex element
     * 
     * @return a pointer to the vertex element
     */    
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;
    
    /**
     * @param pNode pointer to the node
     * 
     * @return a vector of pointers to those vertex elements that are associated with a given node.
     */
    std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> GetElementsOwnedByNode(Node<SPACE_DIM>* pNode);
    
    /**
     * Delete mNodes and mElements.
     */
    void Clear();
    
    /**
     * Re-mesh the mesh.
     * 
     * @param map a NodeMap which associates the indices of VertexElements in the old mesh
     *            with indices of VertexElements in the new mesh.  This should be created 
     *            with the correct size, GetNumElements()
     */
    void ReMesh(NodeMap& elementMap);

    /**
     * Helper method for ReMesh to perform the T-1 Swap
     * @param NodeA one of the nodes to perform the swap with 
     *        NodeB the other node to perform the swap
     */  
    void T1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);
    
//    /**
//     * Alternative version of the ReMesh() method, which does not require a NodeMap. 
//     * 
//     * Note: inherited classes should overload ReMesh(NodeMap&).
//     */
//    void ReMesh();
    
    // when will these be needed?
//    unsigned GetNumAllNodes() const;
//    unsigned GetNumAllElements();

//    
//    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
//                                         bool cullInternalFaces=false)
//    {}
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, 
                                               std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements)
{
    Clear();
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<SPACE_DIM>* temp_node = nodes[index];
        mNodes.push_back(temp_node);
    }
    
    for (unsigned index=0; index<vertexElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* temp_vertex_element = vertexElements[index];
        mElements.push_back(temp_vertex_element);
    }
    
    SetupVertexElementsOwnedByNodes();
    
    mAllocatedMemory = false;
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexMesh(unsigned numAcross, unsigned numUp)
{
    assert(numAcross>1);
    unsigned node_index = 0;
    
    // Create the nodes
    for (unsigned j=0;j<=2*numUp+1;j++)
    {
        if (j%2 == 0)
        {
            for (unsigned i=1;i<=3*numAcross+1;i+=2)
            {
                if (j!=0 || i!= 3*numAcross+1)
                {
                    if (i%3 != 2)
                    {
                        Node<2>* p_node = new Node<2>(node_index, false, i/(2.0*sqrt(3)),j/2.0);
                        mNodes.push_back(p_node);
                        node_index++;
                    }
                }
            }
        }
        else 
        {
            for (unsigned i=0;i<=3*numAcross+1;i+=2)
            {
                if ((j!=2*numUp+1 || i != 0) && (j!=2*numUp+1 || i!= 3*numAcross+1))
                {
                    if (i%3 != 2)
                    {
                        Node<2>* p_node = new Node<2>(node_index, false, i/(2.0*sqrt(3)),j/2.0);
                        mNodes.push_back(p_node);
                        node_index++;
                    }
                }
            }
        }
    }  
    
    // Create the elements. The array node_indices contains the 
    // global node indices from bottom left, going anticlockwise.
    
    unsigned node_indices[6];
    unsigned element_index;
    
    for (unsigned j=0;j<numUp;j++)
    {
        {
            for (unsigned i=0; i<numAcross;i++)
            {
                element_index=j*numAcross+i;
                
                if (numAcross%2==0) // numAcross is even
                {
                    if (j == 0)     // Bottom row
                    {
                        if (i%2 == 0) // even
                        {
                            node_indices[0] = i;
                        }
                        else // odd
                        {
                            node_indices[0] = numAcross+i;
                        }                                           
                    }                       
                    else    // not on the bottom row 
                    {
                         if (i%2 == 0) // even
                        {
                            node_indices[0] = (2*numAcross+1)+2*(j-1)*(numAcross+1)+i;
                        }
                        else // odd
                        {
                            node_indices[0] = (2*numAcross+1)+(2*j-1)*(numAcross+1)+i;
                        }                        
                    }
                        
                }
                else // numAcross is odd
                {
                    if (i%2 == 0) //Even
                    {
                        node_indices[0] = 2*j*(numAcross+1)+i;
                    }
                    else // odd
                    {
                        node_indices[0] = (2*j+1)*(numAcross+1)+i;
                    }
                }
                node_indices[1] = node_indices[0]+1;
                node_indices[2] = node_indices[0]+numAcross+2;
                node_indices[3] = node_indices[0]+2*numAcross+3;
                node_indices[4] = node_indices[0]+2*numAcross+2;
                node_indices[5] = node_indices[0]+numAcross+1;
                 
                if ((j==numUp-1)&&(i%2 == 1))
                {
                    // On top row and its an odd column nodes 
                    node_indices[3]-=1;
                    node_indices[4]-=1;
                }
                  
                if ((j==0)&&(i%2 == 0)&&(numAcross%2==0))
                {
                    // On bottom row and its an even column and there is
                    // an even number of columns in total, (i.e. the very bottom) 
                    node_indices[2]-=1;
                    node_indices[3]-=1;
                    node_indices[4]-=1;
                    node_indices[5]-=1;
                }

                std::vector<Node<2>*> element_nodes;
                
                for (int i=0; i<6; i++)
                {
                   element_nodes.push_back(mNodes[node_indices[i]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                mElements.push_back(p_element);
            }
        }
    }  

//    HoneycombMeshGenerator generator(numAcross+1,numUp+1,0,false);
//    MutableMesh<2,2>* p_mesh = generator.GetMesh();
//    VoronoiTessellation<2> tessellation(*p_mesh);
//    
//    for (unsigned i = 0;i<tessellation.GetNumVertices();i++)
//    {
//        c_vector<double,2>* position = tessellation.GetVertex(i);
//        Node<2>* p_node = new Node<2>(0, false, (*position)(0), (*position)(1));
//        mNodes.push_back(p_node);
//    }    

    /// \todo: loop over the p_mesh's nodes, and if it is a non-boundary node create a VertexElement using
    //         the corresponding cell. Then get rid of the nodes in mNodes that do not belong in any cell.

    mAllocatedMemory = true;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::~VertexMesh()
{
    if (mAllocatedMemory)
    {
        Clear();
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::SetupVertexElementsOwnedByNodes()
{
    for (unsigned index=0; index<mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }
};


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }

    mNodes.clear();
    mElements.clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNodes.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index) const
{
    assert(index < mNodes.size());
    return mNodes[index];
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM,SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(NodeMap& elementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 || SPACE_DIM==3 );
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE

    // Make sure the map is big enough
    elementMap.Resize(GetNumElements());
    
    if (SPACE_DIM==2)
    {
        /// \todo put code for remeshing in 2D here (see #827)
    
        unsigned new_index = 0;
        for (unsigned i=0; i<GetNumElements(); i++)
        {
            if (mElements[i]->IsDeleted())
            {
                elementMap.SetDeleted(i);
            }
            else
            {
                elementMap.SetNewIndex(i, new_index);
                new_index++;
            }
        }
        
        /*
         * We do not need to call Clear() and remove all current data, since
         * cell birth, rearrangement and death result only in local remeshing
         * of a vertex-based mesh.
         * 
         * Instead, we should now remove any deleted nodes and elements.
         * 
         * We should then construct any new nodes, including boundary nodes; 
         * then new elements; then new edges.
         * 
         * Finally (or should this be at the start?), we should perform any 
         * cell rearrangements.
         */
    }
    else // 3D
    {
        /// \todo put code for remeshing in 3D here (see #827)
    }    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMesh<ELEMENT_DIM, SPACE_DIM>::T1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert( SPACE_DIM==2 ); // only works in 2D at present
    assert( ELEMENT_DIM == SPACE_DIM );
    #undef COVERAGE_IGNORE

    // Find elements containing nodes A and B
    
    // Move Nodes A to C and node B to D.
    
    // Restructure elements -- Remember to update nodes and elements
}

//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllNodes()
//{
//    return mNodes.size();
//}
//
//
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//unsigned VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements()
//{
//    return mElements.size();
//}
//
   

#endif /*VERTEXMESH_HPP_*/
