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

#ifndef ABSTRACTMESH_HPP_
#define ABSTRACTMESH_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>

#include "Node.hpp"
#include "BoundaryElement.hpp"
#include "Element.hpp"
#include "AbstractMeshReader.hpp"
#include "DistributedVectorFactory.hpp"


/**
 * Abstract base class for all meshes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractMesh
{
private:

    /**
     * Pure virtual solve node mapping method.
     * Overridden in TetrahedralMesh and ParallelTetrahedralMesh classes.
     *
     * @param index
     */
    virtual unsigned SolveNodeMapping(unsigned index) const = 0;

    /**
     * Pure virtual solve element mapping method.
     * Overridden in TetrahedralMesh and ParallelTetrahedralMesh classes.
     *
     * @param index
     */
    virtual unsigned SolveElementMapping(unsigned index) const = 0;

    /**
     * Pure virtual solve boundary element mapping method.
     * Overridden in TetrahedralMesh and ParallelTetrahedralMesh classes.
     *
     * @param index
     */
    virtual unsigned SolveBoundaryElementMapping(unsigned index) const = 0;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the mesh.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       // Don't do anything - this is just so subclasses can archive member variables.
    }

protected:  // Give access of these variables to subclasses

    /** Vector of pointers to nodes in the mesh. */
    std::vector<Node<SPACE_DIM> *> mNodes;

    /** Vector of pointers to boundary nodes in the mesh. */
    std::vector<Node<SPACE_DIM> *> mBoundaryNodes;

    /** Vector of pointers to elements in the mesh. */
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> mElements;

    /** Vector of pointers to boundary elements in the mesh. */
    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;

    /** Vector containing the number of nodes owned by each processor.  Redundant #988*/
    std::vector<unsigned> mNodesPerProcessor;

    /** DistributedVectorFactory capable of reproducing the given number of nodes owned by each processor. */
    DistributedVectorFactory *mpDistributedVectorFactory;

    /** Vector containing node permutation information. */
    std::vector<unsigned> mNodesPermutation;

    /**
     * If the mesh is constructed from file using a MeshReader, this member
     * variable stores the base name of these files.
     */
    std::string mMeshFileBaseName;

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Definition of element Iterator type. */
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::const_iterator ElementIterator;
    /** Definition of boundary element Iterator type. */
    typedef typename std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *>::const_iterator BoundaryElementIterator;
    /** Definition of boundary node Iterator type. */
    typedef typename std::vector<Node<SPACE_DIM> *>::const_iterator BoundaryNodeIterator;

    /** Forward declaration */
    class NodeIterator;
    
    /**
     * Get an iterator to the first node in the mesh.
     * 
     * @param skipDeletedNodes whether to include deleted nodes
     */
    inline NodeIterator GetNodeIteratorBegin(bool skipDeletedNodes=true);
    
    /**
     * Get an iterator to one past the last node in the mesh.
     */
    inline NodeIterator GetNodeIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Constructor.
     */
    AbstractMesh();

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractMesh();

    /**
     * Get the number of nodes that are actually in use.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * Get the number of elements that are actually in use.
     */
    virtual unsigned GetNumElements() const;

    /**
     * Get the number of boundary elements that are actually in use.
     */
    virtual unsigned GetNumBoundaryElements() const;

    /**
     * Get the number of boundary nodes in the mesh.
     */
    unsigned GetNumBoundaryNodes();// should this be overloaded and virtual too?

    /**
     * Get the total number of nodes (including those marked as deleted).
     */
    unsigned GetNumAllNodes() const;

    /**
     * Get the total number of elements (including those marked as deleted).
     */
    unsigned GetNumAllElements();

    /**
     * Get the total number of boundary elements (including those marked as deleted).
     */
    unsigned GetNumAllBoundaryElements();

    /**
     * Get the node with a given index in the mesh.
     *
     * @param index
     * @return a pointer to the node.
     */
    Node<SPACE_DIM> *GetNode(unsigned index) const;

    /**
     * Get the element with a given index in the mesh.
     *
     * @param index
     * @return a pointer to the element.
     */
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Get the boundary element with a given index in the mesh.
     *
     * @param index
     * @return a pointer to the boundary element.
     */
    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* GetBoundaryElement(unsigned index) const;

    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     *
     * @param lo is the lowest node number owned by the process
     * @param hi is one higher than the highest node number owned by the process
     * ie. this process owns nodes [lo..hi)
     * and element is "owned" if one or more of its nodes are owned
     */
    virtual void SetElementOwnerships(unsigned lo, unsigned hi);

    /**
     * Construct the mesh using a MeshReader.
     * This method must be overridden in concrete classes.
     *
     * @param rMeshReader the mesh reader
     * @param cullInternalFaces whether to cull internal faces (defaults to false)
     */
    virtual void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader,
                                         bool cullInternalFaces=false)=0;

    /**
     * Read in the number of nodes per processor from file.
     *
     * @param nodesPerProcessorFile
     */
    virtual void ReadNodesPerProcessorFile(const std::string& nodesPerProcessorFile);

//Going
    /**
     * Get method for mNodesPerProcessor.
     */
    std::vector<unsigned>& rGetNodesPerProcessor();

    /**
     * Get method for mNodesPerProcessor.
     */
    DistributedVectorFactory * GetDistributedVectorFactory();

    /**
     * Permute the nodes so that they appear in a different order in mNodes
     * (and their mIndex's are altered accordingly).
     */
    virtual void PermuteNodes();

    /**
     * Return a pointer to the first element in the mesh.
     */
    ElementIterator GetElementIteratorBegin() const;

    /**
     * Return a pointer to *one past* the last element in the mesh
     * (for consistency with STL iterators).
     */
    ElementIterator GetElementIteratorEnd() const;

    /**
     * Return a pointer to the first boundary element in the mesh.
     */
    BoundaryElementIterator GetBoundaryElementIteratorBegin() const;

    /**
     * Return a pointer to *one past* the last boundary element in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryElementIterator GetBoundaryElementIteratorEnd() const;

    /**
     * Return a pointer to the first boundary node in the mesh.
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorBegin() const;

    /**
     * Return a pointer to *one past* the last boundary node in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorEnd() const;

    /**
     * Compute the inverse Jacobian for a given element in the mesh.
     *
     * @param elementIndex index of an element
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian matrix
     * @param rInverseJacobian  the inverse Jacobian matrix
     */
    virtual void GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double &rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const;

    /**
     * Compute the weighted direction for a given boundary element.
     *
     * @param elementIndex index of an element
     * @param rWeightedDirection the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian matrix
     */
    virtual void GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double &rJacobianDeterminant) const;

    /**
     * Get method for mMeshFileBaseName.
     */
    std::string GetMeshFileBaseName() const;

    /**
     * Get method for mNodesPermutation.
     */
    std::vector<unsigned>& rGetNodePermutation();
    
    //////////////////////////////////////////////////////////////////////
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////
    
    /**
     * A smart iterator over the nodes in the mesh.
     */
    class NodeIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current node.
         * 
         * Make sure to use a reference for the result to avoid copying nodes unnecessarily.
         */
        inline Node<SPACE_DIM>& operator*();
        
        /**
         * Member access from a pointer.
         */
        inline Node<SPACE_DIM>* operator->();
        
        /**
         * Comparison not-equal-to.
         *
         * @param other iterator with which comparison is made
         */
        inline bool operator!=(const AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator& other);
        
        /**
         * Prefix increment operator.
         */
        inline NodeIterator& operator++();
        
        /**
         * Constructor for a new iterator.
         * 
         * This should not be called directly by user code; use the mesh methods
         * AbstractMesh::GetNodeIteratorBegin and AbstractMesh::GetNodeIteratorEnd instead.
         * 
         * @param rMesh the mesh to iterator over
         * @param nodeIter where to start iterating
         * @param skipDeletedNodes whether to include deleted nodes
         */
        NodeIterator(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                     typename std::vector<Node<SPACE_DIM> *>::iterator nodeIter,
                     bool skipDeletedNodes=true);
    private:
        /** The mesh we're iterating over. */
        AbstractMesh& mrMesh;

        /** The actual node iterator. */
        typename std::vector<Node<SPACE_DIM> *>::iterator mNodeIter;
        
        /** Whether to skip deleted nodes. */
        bool mSkipDeletedNodes;
        
        /**
         * Helper method to say when we're at the end.
         */
        inline bool IsAtEnd();
        
        /**
         * Helper method to say if we're allowed to point at this node.
         */
        inline bool IsAllowedNode();
    };
};

namespace boost
{
namespace serialization
{
/**
 * Since this abstract class is templated, we cannot use
 * the preprocessor macro BOOST_IS_ABSTRACT, and instead
 * must drop down to the underlying source code.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct is_abstract<AbstractMesh<ELEMENT_DIM, SPACE_DIM> >
{
    /** The type that is an abstract class. */
    typedef mpl::bool_<true> type;
    /** The type is an abstract class, so value=true. */
    BOOST_STATIC_CONSTANT(bool, value=true);
};
}
}

//////////////////////////////////////////////////////////////////////////////
//         Iterator class implementation - most methods are inlined         //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeIteratorBegin(
        bool skipDeletedNodes)
{
    return NodeIterator(*this, mNodes.begin(), skipDeletedNodes);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeIteratorEnd()
{
    return NodeIterator(*this, mNodes.end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator*()
{
    assert(!IsAtEnd());
    return **mNodeIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator->()
{
    assert(!IsAtEnd());
    return *mNodeIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator!=(const AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator& other)
{
    return mNodeIter != other.mNodeIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator& AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::operator++()
{
    do
    {
        ++mNodeIter;
    }
    while (!IsAtEnd() && !IsAllowedNode());
    
    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::NodeIterator(
        AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        typename std::vector<Node<SPACE_DIM> *>::iterator nodeIter,
        bool skipDeletedNodes)
    : mrMesh(rMesh),
      mNodeIter(nodeIter),
      mSkipDeletedNodes(skipDeletedNodes)
{
    if (mrMesh.mNodes.size() == 0)
    {
        // Cope with empty meshes
        mNodeIter = mrMesh.mNodes.end();
    }
    else
    {
        // Make sure we start at an allowed node
        if (mNodeIter == mrMesh.mNodes.begin() && !IsAllowedNode())
        {
            ++(*this);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::IsAtEnd()
{
    return mNodeIter == mrMesh.mNodes.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator::IsAllowedNode()
{
    return !(mSkipDeletedNodes && (*this)->IsDeleted());
}



#endif /*ABSTRACTMESH_HPP_*/
