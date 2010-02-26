/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef MESHBASEDTISSUEWITHGHOSTNODES_HPP_
#define MESHBASEDTISSUEWITHGHOSTNODES_HPP_

#include "MeshBasedTissue.hpp"
#include "TrianglesMeshReader.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A facade class encapsulating a mesh-based 'tissue' with ghost nodes.
 *
 * If simulating a crypt with a mesh-based tissue, the mesh should be surrounded by at
 * least one layer of ghost nodes. These are nodes which do not correspond to a cell,
 * but are necessary for remeshing (because the remesher tries to create a convex hull
 * of the set of nodes) and visualization purposes. The MeshBasedTissueWithGhostNodes
 * class deals with these ghost nodes, hiding the 'ghost nodes' concept from the
 * TissueSimulation class, so the latter only ever deals with real cells.
 */
template<unsigned DIM>
class MeshBasedTissueWithGhostNodes : public MeshBasedTissue<DIM>
{
private:
    /** Just so that the test can test the private functions */
    friend class TestMeshBasedTissueWithGhostNodes;

    /** Records whether a node is a ghost node or not */
    std::vector<bool> mIsGhostNode;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This needs to be first so that MeshBasedTissue::Validate() doesn't go mental.
        archive & mIsGhostNode;
        archive & boost::serialization::base_object<MeshBasedTissue<DIM> >(*this);
    }

    /**
     * Set the ghost nodes by taking in a set of which nodes indices are ghost nodes.
     *
     * @param rGhostNodeIndices set of node indices corresponding to ghost nodes
     */
    void SetGhostNodes(const std::set<unsigned>& rGhostNodeIndices);

    /**
     * This is called after a tissue has been constructed to check the
     * user gave consistent instructions. Check consistency of our
     * internal data structures:
     * Each node must have a cell associated with it OR must be a ghost node.
     *
     * It is called after cells are added or removed from MeshBasedTissue
     * as it is an overridden virtual method.
     */
    void Validate();

public:

    /**
     * Create a new tissue facade from a mesh and collection of cells.
     *
     * @param rMesh a mutable tetrahedral mesh
     * @param rCells TissueCells corresponding to the nodes of the mesh
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the tissue to free the mesh memory on destruction
     */
    MeshBasedTissueWithGhostNodes(MutableMesh<DIM, DIM>& rMesh,
                                  const std::vector<TissueCell>& rCells,
                                  const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                                  bool deleteMesh=false);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable tetrahedral mesh.
     */
    MeshBasedTissueWithGhostNodes(MutableMesh<DIM, DIM>& rMesh);

    /**
     * Overridden UpdateNodeLocation() method.
     *
     * Update the location of each node in the tissue given
     * a vector of forces on nodes and a time step over which
     * to integrate the equations of motion.
     *
     * @param rNodeForces  forces on nodes
     * @param dt  time step
     */
    void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt);

    /**
     * @return mIsGhostNode.
     */
    std::vector<bool>& rGetGhostNodes();

    /**
     * Overridden IsGhostNode() method.
     *
     * Find if a given node is a ghost node. The abstract method always returns false
     * but is overridden in subclasses.
     *
     * @param index the global index of a specified node
     *
     * @return whether the node is a ghost node
     */
    bool IsGhostNode(unsigned index);

    /**
     * @return the indices of those nodes that are ghost nodes.
     */
    std::set<unsigned> GetGhostNodeIndices();

    /**
     * Update the GhostNode positions using the spring force model with rest length=1.
     * Forces are applied to ghost nodes from connected ghost and normal nodes.
     *
     * @param dt
     */
    void UpdateGhostPositions(double dt);

    /**
     * Update mIsGhostNode if required by a remesh.
     *
     * @param rMap A map between node indices before and after remesh
     */
    void UpdateGhostNodesAfterReMesh(NodeMap& rMap);

    /**
     * This method is used to calculate the force between GHOST nodes.
     *
     * @param rNodeAGlobalIndex
     * @param rNodeBGlobalIndex
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex);

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the tissue and update mIsGhostNode.
     *
     * @param rNewCell  the cell to add
     * @param rCellDivisionVector  the position in space at which to put it
     * @param pParentCell pointer to a parent cell (if required)
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell* AddCell(TissueCell& rNewCell, const c_vector<double,DIM>& rCellDivisionVector, TissueCell* pParentCell=NULL);

    /**
     * Overridden GenerateCellResults() method.
     *  Generate results for a given cell in the current tissue state to output files.
     *
     * @param locationIndex location index of the cell
     * @param rCellProliferativeTypeCounter cell type counter
     * @param rCellMutationStateCounter cell mutation state counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void GenerateCellResults(unsigned locationIndex,
                             std::vector<unsigned>& rCellProliferativeTypeCounter,
                             std::vector<unsigned>& rCellMutationStateCounter,
                             std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Overridden GenerateCellResultsAndWriteToFiles() method.
     *
     * Call GenerateCellResults() on each cell then call WriteCellResultsToFiles().
     * Also accounts for ghost nodes.
     */
    virtual void GenerateCellResultsAndWriteToFiles();
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedTissueWithGhostNodes)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MeshBasedTissueWithGhostNodes.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MeshBasedTissueWithGhostNodes<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const MutableMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a MeshBasedTissueWithGhostNodes.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MeshBasedTissueWithGhostNodes<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    MutableMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)MeshBasedTissueWithGhostNodes<DIM>(*p_mesh);

}
}
} // namespace

#endif /*MESHBASEDTISSUEWITHGHOSTNODES_HPP_*/
