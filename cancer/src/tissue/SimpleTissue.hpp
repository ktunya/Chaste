#ifndef SIMPLETISSUE_HPP_
#define SIMPLETISSUE_HPP_

#include "AbstractTissue.cpp"
#include "ConformingTetrahedralMesh.cpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

template<unsigned DIM>
class SimpleTissue : public AbstractTissue<DIM>
{
    friend class TestSimpleTissue;
private:

    /** List of nodes */
    std::vector<Node<DIM> > mNodes;
    
    
    friend class boost::serialization::access;
    /**
     * Serialize the facade.
     * 
     * Note that serialization of the nodes is handled by load/save_construct_data,
     * so we don't actually have to do anything here except delegate to the base class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTissue<DIM> >(*this);
        
        Validate(); // paranoia
    }
    
            
    /** 
     * Add a new node to the tissue. 
     */
    unsigned AddNode(Node<DIM> *pNewNode);
    
    /** 
     * Update the correspondence between nodes and cells.
     */
    void UpdateNodeCellMap();

    /** 
     * Remove the node with a given index.
     */
    void RemoveNode(unsigned index);
    
    /** 
     * Move the node with a given index to a new point in space.
     */
    void SetNode(unsigned index, ChastePoint<DIM> point);
    

public:
    SimpleTissue(const std::vector<Node<DIM> >& rNodes, const std::vector<TissueCell>& rCells);

    /**
     * Constructor for use by the archiving - doesn't take in cells, since these are
     * dealt with by the serialize method of our base class.
     */
    SimpleTissue(const std::vector<Node<DIM> >& rNodes);
    
    ~SimpleTissue() 
    {}

    /** 
     * Get the number of nodes in the tissue.
     */
    unsigned GetNumNodes();
    
    /**
     * Get a pointer to the node with a given index.
     */  
    Node<DIM>* GetNode(unsigned index);
    
    /**
     * Get a pointer to the node corresponding to a given TissueCell.
     */
    Node<DIM>* GetNodeCorrespondingToCell(const TissueCell& rCell);

    /**
     * Add a new cell to the tissue.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell* AddCell(TissueCell cell, c_vector<double,DIM> newLocation);
    
    /** 
     * Remove all cells labelled as dead. 
     * 
     * Note that after calling this method the tissue will be in an inconsistent state until 
     * the equivalent of a 'remesh' is performed! So don't try iterating over cells or anything 
     * like that.
     * 
     *  @return number of cells removed
     */
    unsigned RemoveDeadCells();
    
    /**
     * Check consistency of our internal data structures.
     */
    void Validate();
    
    void WriteResultsToFiles(bool OutputCellTypes);
    
    std::vector<Node<DIM> >& rGetNodes();
    const std::vector<Node<DIM> >& rGetNodes() const;
    
    /**
     * Move a cell to a new location.
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    void MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation);
    
};

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleTissue)

namespace boost {
namespace serialization {

/**
 * Non-intrusive serialization for Node - save method.
 */
template<class Archive, unsigned SPACE_DIM>
inline void save(
    Archive & ar,
    const Node<SPACE_DIM> &rNode,
    const unsigned int /* file_version */)
{
    // Save deleted flag
    const bool is_deleted = rNode.IsDeleted();
    ar << is_deleted;
}

/**
 * Non-intrusive serialization for Node - load method.
 */
template<class Archive, unsigned SPACE_DIM>
inline void load(
    Archive & ar,
    Node<SPACE_DIM> &rNode,
    const unsigned int /* file_version */)
{
    // Load deleted flag
    bool is_deleted;
    ar >> is_deleted;
    #define COVERAGE_IGNORE
    if (is_deleted)
    {
        rNode.MarkAsDeleted();
    }
    #undef COVERAGE_IGNORE
}


/**
 * Non-intrusive serialization for Node - serialize method.
 * This calls save or load as appropriate.
 */
template<class Archive, unsigned SPACE_DIM>
inline void serialize(
    Archive & ar,
    Node<SPACE_DIM>& rNode,
    const unsigned int file_version)
{
    boost::serialization::split_free(ar, rNode, file_version);
}


/**
 * Serialize information required to construct a Node.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const Node<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    const unsigned index = t->GetIndex();
    ar << index;
    const bool is_boundary = t->IsBoundaryNode();
    ar << is_boundary;
    const c_vector<double, DIM>& r_loc = t->rGetLocation();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << r_loc[i];
    }
}

/**
 * De-serialize constructor parameters and initialise a Node.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, Node<DIM> * t, const unsigned int file_version)
{
    unsigned index;
    ar >> index;
    bool is_boundary;
    ar >> is_boundary;
    // Load the location
    c_vector<double, DIM> loc;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> loc[i];
    }
    
    // Invoke inplace constructor to initialize instance
    ::new(t)Node<DIM>(index, loc, is_boundary);
}


/**
 * Serialize information required to construct a SimpleTissue facade.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SimpleTissue<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    //const ConformingTetrahedralMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & t->rGetNodes();
}

/**
 * De-serialize constructor parameters and initialise SimpleTissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SimpleTissue<DIM> * t, const unsigned int file_version)
{
    // Load the nodes
    std::vector<Node<DIM> > nodes;
    ar >> nodes;
    
    // Invoke inplace constructor to initialize instance
    ::new(t)SimpleTissue<DIM>(nodes);
}

}} // close namespaces


#endif /*SIMPLETISSUE_HPP_*/
