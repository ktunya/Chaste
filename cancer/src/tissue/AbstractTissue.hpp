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
#ifndef ABSTRACTTISSUE_HPP_
#define ABSTRACTTISSUE_HPP_

#include "TissueCell.hpp"
#include "OutputFileHandler.hpp"

#include <list>

#include <climits> // work around boost bug
#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>


/**
 * An abstract facade class encapsulating a tissue.
 *
 * Contains a group of cells and associated methods.
 *
 */
template<unsigned DIM>
class AbstractTissue
{
protected:

    /** List of cells */
    std::list<TissueCell> mCells;

    /** Map location (node or VertexElement) indices back to cells */
    std::map<unsigned, TissueCell*> mLocationCellMap;

    /** Current cell mutation state counts */
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> mCellMutationStateCount;

    /** Current cell type counts */
    c_vector<unsigned, NUM_CELL_TYPES> mCellTypeCount;

    /** Current cell cycle phase counts */
    c_vector<unsigned, 5> mCellCyclePhaseCount;

    /** Results file for node visualization */
    out_stream mpVizNodesFile;
    
    /** Results file for cell visualization */
    out_stream mpVizCellTypesFile;

    /** Results file for cell mutation states */
    out_stream mpCellMutationStatesFile;

    /** Results file for cell ancestors */
    out_stream mpCellAncestorsFile;

    /** Results file for cell types */
    out_stream mpCellTypesFile;

    /** Results file for cell cycle phases */
    out_stream mpCellCyclePhasesFile;

    /** Results file for cell variables */
    out_stream mpCellVariablesFile;

    /** Whether the tissue contains a mesh */
    bool mTissueContainsMesh;

    /** Whether the tissue contains ghost nodes */
    bool mTissueContainsGhostNodes;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCells;
        archive & mLocationCellMap;
        archive & mTissueContainsMesh;
        archive & mTissueContainsGhostNodes;
    }

public:

    /**
     * Default constructor.
     */
    AbstractTissue(const std::vector<TissueCell>& rCells);

    /**
     * Constructor for use by archiving - doesn't take in cells, since these are dealt
     * with by the serialize method.
     */
    AbstractTissue()
    {}

    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractTissue()
    {}

    /**
     * Initialise each cell's cell cycle model.
     */
    void InitialiseCells();

    /**
     * @return reference to mCells.
     */
    std::list<TissueCell>& rGetCells();
    
    /**
     * @return const reference to mCells (used in archiving).
     */
    const std::list<TissueCell>& rGetCells() const;

    /**
     * @return whether the tissue contains a mesh.
     */
    bool HasMesh();

    /**
     * @return whether the tissue contains ghost nodes.
     */
    bool HasGhostNodes();

    /**
     * @return the number of nodes in the tissue.
     * 
     * As this method is pure virtual, it must be overridden 
     * in subclasses.
     */
    virtual unsigned GetNumNodes()=0;

    /**
     * @return a pointer to the node with a given index.
     * 
     * As this method is pure virtual, it must be overridden 
     * in subclasses.
     * 
     * @param index  global index of the specified node
     */
    virtual Node<DIM>* GetNode(unsigned index)=0;
    
    /**
     * Add a new node to the tissue.
     * 
     * As this method is pure virtual, it must be overridden 
     * in subclasses.
     * 
     * @param pNewNode pointer to the new node 
     * @return global index of new node in tissue
     */
    virtual unsigned AddNode(Node<DIM> *pNewNode)=0;

    /**
     * Move the node with a given index to a new point in space.
     * 
     * As this method is pure virtual, it must be overridden 
     * in subclasses.
     * 
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    virtual void SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)=0;
    
    /**
     * @return whether a given cell is associated with a ghost node.
     */
    virtual bool IsCellAssociatedWithAGhostNode(TissueCell cell);
    
    /**
     * @return whether a given cell is associated with a deleted node.
     * 
     * As this method is pure virtual, it must be overridden 
     * in subclasses.
     */
    virtual bool IsCellAssociatedWithADeletedNode(TissueCell cell)=0;

    /**
     * Get the damping constant for this node - ie d in drdt = F/d.
     * 
     * As this method is pure virtual, it must be overridden 
     * in subclasses.
     * 
     * @param nodeIndex the global index of this node
     * @return the damping constant at the node
     */
    virtual double GetDampingConstant(unsigned nodeIndex)=0;

    /**
     * Add a new cell to the tissue.
     * 
     * As this method is pure virtual, it must be overridden 
     * in subclasses.
     *
     * @param rNewCell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @param pParentCell pointer to a parent cell (if required)
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual TissueCell* AddCell(TissueCell& rNewCell, c_vector<double,DIM> newLocation, TissueCell* pParentCell=NULL)=0; 

    class Iterator; // Forward declaration; see below

    /**
     * Remove all cells labelled as dead.
     * 
     * As this method is pure virtual, it must be overridden 
     * in subclasses.
     * 
     * @return number of cells removed
     */
    virtual unsigned RemoveDeadCells()=0;

    /**
     * Remove the Nodes (for cell-centre) or VertexElements (for cell-vertex) which 
     * have been marked as deleted and update the correspondence with TissueCells.
     */
    virtual void Update()=0;

    /**
     * Check consistency of our internal data structures.
     */
    virtual void Validate()=0;

    /**
     * Find out how many cells of each mutation state there are
     *
     * @return The number of cells of each mutation state (evaluated at each visualizer output)
     * [0] = healthy count
     * [1] = labelled cells
     * [2] = APC one hit
     * [3] = APC two hit
     * [4] = beta catenin one hit
     */
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> GetCellMutationStateCount();

     /**
     * Find out how many cells of each type there are
     *
     * @return The number of cells of each type (evaluated at each visualizer output)
     * [0] = STEM
     * [1] = TRANSIT
     * [2] = DIFFERENTIATED
     * [3] = APOPTOTIC
     */
    c_vector<unsigned, NUM_CELL_TYPES> GetCellTypeCount();

    /**
     * Find out how many cells in each cell cycle phase there are
     *
     * @return The number of cells of each phase (evaluated at each visualizer output)
     * [0] = G_ZERO_PHASE
     * [1] = G_ONE_PHASE
     * [2] = S_PHASE
     * [3] = G_TWO_PHASE
     * [4] = M_PHASE
     */
    c_vector<unsigned, 5> GetCellCyclePhaseCount();

    /**
     * Find if a given node is a ghost node. The abstract method always returns false 
     * but is overridden in subclasses.
     * 
     * @param index the global index of a specified node
     * 
     * @return whether the node is a ghost node
     */
    virtual bool IsGhostNode(unsigned index);

    /**
     * Get the number of real cells.
     */
    unsigned GetNumRealCells();

    /**
     * Sets the Ancestor index of all the cells at this time to be the
     * same as their node index, can be used to trace clonal populations.
     */
    void SetCellAncestorsToNodeIndices();

    /**
     * Loops over cells and makes a list of the ancestors that
     * are part of the tissue.
     *
     * @return remaining_ancestors  The size of this set tells you how many clonal populations remain.
     */
    std::set<unsigned> GetCellAncestors();

    /**
     *  Get the cell corresponding to a given node.
     *
     *  Currently assumes there is one cell for each node, and they are ordered identically in their vectors.
     *  An assertion fails if not.
     */
    TissueCell& rGetCellUsingLocationIndex(unsigned index);
    
    /**
     * If the tissue contains a mesh, write this to file. For use by 
     * the TissueSimulationArchiver. Must be overridden in each subclass
     * that contains a mesh.
     * 
     * @param rArchiveDirectory directory in which archive is stored
     * @param rMeshFileName base name for mesh files
     */
    virtual void WriteMeshToFile(const std::string &rArchiveDirectory, const std::string &rMeshFileName);
        

    /**
     * Use an output file handler to create output files for visualizer and post-processing.
     * 
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param rCleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     * @param outputCellMutationStates  whether to create a cell mutation state results file
     * @param outputCellTypes  whether to create a cell type results file
     * @param outputCellVariables  whether to create a cell-cycle variable results file
     * @param outputCellCyclePhases  whether to create a cell-cycle phase results file
     * @param outputCellAncestors  whether to create a cell ancestor results file
     */
    virtual void CreateOutputFiles(const std::string &rDirectory,
                                   bool rCleanOutputDirectory,
                                   bool outputCellMutationStates,
                                   bool outputCellTypes,
                                   bool outputCellVariables,
                                   bool outputCellCyclePhases,
                                   bool outputCellAncestors);

    virtual void WriteResultsToFiles(bool outputCellMutationStates,
                                     bool outputCellTypes,
                                     bool outputCellVariables,
                                     bool outputCellCyclePhases,
                                     bool outputCellAncestors)
    {}

    /**
     * Write the current time and node results to output files.
     * 
     * @param outputCellMutationStates  whether to output cell mutation state results
     * @param outputCellTypes  whether to output cell type results
     * @param outputCellVariables  whether to output cell-cycle variable results
     * @param outputCellCyclePhases  whether to output cell-cycle phase results
     * @param outputCellAncestors  whether to output cell ancestor results
     * @param rCellTypeCounter cell type counter
     * @param rCellMutationStateCounter cell mutation state counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void WriteTimeAndNodeResultsToFiles(bool outputCellMutationStates,
                                        bool outputCellTypes,
                                        bool outputCellVariables,
                                        bool outputCellCyclePhases,
                                        bool outputCellAncestors,
                                        std::vector<unsigned>& rCellTypeCounter,
                                        std::vector<unsigned>& rCellMutationStateCounter,
                                        std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Generate results for a given cell in the current tissue state to output files.
     * 
     * @param locationIndex location index of the cell
     * @param outputCellMutationStates  whether to output cell mutation state results
     * @param outputCellTypes  whether to output cell type results
     * @param outputCellVariables  whether to output cell-cycle variable results
     * @param outputCellCyclePhases  whether to output cell-cycle phase results
     * @param outputCellAncestors  whether to output cell ancestor results
     * @param rCellTypeCounter cell type counter
     * @param rCellMutationStateCounter cell mutation state counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void GenerateCellResults(unsigned locationIndex,
                             bool outputCellMutationStates,
                             bool outputCellTypes,
                             bool outputCellVariables,
                             bool outputCellCyclePhases,
                             bool outputCellAncestors,
                             std::vector<unsigned>& rCellTypeCounter,
                             std::vector<unsigned>& rCellMutationStateCounter,
                             std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Write the current state of each cell to output files.
     * 
     * @param outputCellMutationStates  whether to output cell mutation state results
     * @param outputCellTypes  whether to output cell type results
     * @param outputCellVariables  whether to output cell-cycle variable results
     * @param outputCellCyclePhases  whether to output cell-cycle phase results
     * @param outputCellAncestors  whether to output cell ancestor results
     * @param rCellTypeCounter cell type counter
     * @param rCellMutationStateCounter cell mutation state counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void WriteCellResultsToFiles(bool outputCellMutationStates,
                                 bool outputCellTypes,
                                 bool outputCellVariables,
                                 bool outputCellCyclePhases,
                                 bool outputCellAncestors,
                                 std::vector<unsigned>& rCellTypeCounter,
                                 std::vector<unsigned>& rCellMutationStateCounter,
                                 std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Close any output files.
     * 
     * @param outputCellMutationStates  whether a cell mutation state results file is open
     * @param outputCellTypes  whether a cell type results file is open
     * @param outputCellVariables  whether a cell-cycle variable results file is open
     * @param outputCellCyclePhases  whether a cell-cycle phase results file is open
     * @param outputCellAncestors  whether a cell ancestor results file is open
     */
    virtual void CloseOutputFiles(bool outputCellMutationStates,
                                  bool outputCellTypes,
                                  bool outputCellVariables,
                                  bool outputCellCyclePhases,
                                  bool outputCellAncestors);

    /**
     * Iterator class allows one to iterate over cells in the tissue.
     * Dereferencing the iterator will give you the current cell.
     * There are also methods to get the node representing this cell,
     * and the location of that node.
     */
    class Iterator
    {
    public:

        /**
         * Dereference the iterator giving you a *reference* to the current cell.
         * Make sure to use a reference for the result to avoid copying cells unnecessarily.
         */
        inline TissueCell& operator*();

        /**
         * Member access from a pointer.
         */
        inline TissueCell* operator->();

        /**
         * Comparison not-equal-to.
         */
        inline bool operator!=(const Iterator& other);

        /**
         * Prefix increment operator.
         */
        inline Iterator& operator++();

        /**
         * Constructor for a new iterator.
         */
        Iterator(AbstractTissue& rTissue, std::list<TissueCell>::iterator cellIter);

        /**
         * The iterator must have a virtual destructor.
         */
        virtual ~Iterator()
        {}

    private:

        /**
         * Private helper function which tells us if we're pointing at a real cell.
         * Assumes we are within range (i.e. not at End).
         *
         * Real cells are not deleted.
         */
        virtual inline bool IsRealCell();

        /**
         * Private helper function saying whether we're at the end of the cells.
         */
        inline bool IsAtEnd();

        /** The tissue member. */
        AbstractTissue& mrTissue;

        /** Cell iterator member. */
        std::list<TissueCell>::iterator mCellIter;

        /** Location index member. */
        unsigned mLocationIndex;
    };

    /**
     * @return iterator pointing to the first cell in the tissue
     */
    Iterator Begin();

    /**
     * @return iterator pointing to one past the last cell in the tissue
     */
    Iterator End();

};

enum cell_colours
{
    STEM_COLOUR, // 0
    TRANSIT_COLOUR, // 1
    DIFFERENTIATED_COLOUR, // 2
    EARLY_CANCER_COLOUR, // 3
    LATE_CANCER_COLOUR, // 4
    LABELLED_COLOUR, // 5
    APOPTOSIS_COLOUR, // 6
    INVISIBLE_COLOUR, // visualizer treats '7' as invisible
};

template<unsigned DIM>
AbstractTissue<DIM>::AbstractTissue(const std::vector<TissueCell>& rCells)
             : mCells(rCells.begin(), rCells.end()),
               mTissueContainsMesh(false),
               mTissueContainsGhostNodes(false)
{
    // There must be at least one cell
    assert(mCells.size() > 0);

    // Set up the map between location indices and cells
    for (std::list<TissueCell>::iterator it = mCells.begin();
         it != mCells.end();
         ++it)
    {
        /// \todo Check it points to a real cell (see #430),
        /// if not do:
        /// it = this->mCells.erase(it); --it; continue; 
        unsigned index = it->GetLocationIndex();
        mLocationCellMap[index] = &(*it);
    }

    // Initialise cell counts to zero
    for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
    {
        mCellMutationStateCount[i] = 0;
    }
    for (unsigned i=0; i<NUM_CELL_TYPES; i++)
    {
        mCellTypeCount[i] = 0;
    }
    for (unsigned i=0; i<5; i++)
    {
        mCellCyclePhaseCount[i] = 0;
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::InitialiseCells()
{
    for (std::list<TissueCell>::iterator iter = mCells.begin();
        iter != mCells.end();
        ++iter)
    {
        iter->InitialiseCellCycleModel();
    }
}

template<unsigned DIM>
std::list<TissueCell>& AbstractTissue<DIM>::rGetCells()
{
    return mCells;
}

template<unsigned DIM>
const std::list<TissueCell>& AbstractTissue<DIM>::rGetCells() const
{
    return this->mCells;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::HasMesh()
{
    return mTissueContainsMesh;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::HasGhostNodes()
{
    return mTissueContainsGhostNodes;
}

template<unsigned DIM>
unsigned AbstractTissue<DIM>::GetNumRealCells()
{
    unsigned counter = 0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        counter++;
    }
    return counter;
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetCellAncestorsToNodeIndices()
{
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->SetAncestor(cell_iter->GetLocationIndex());
    }
}

template<unsigned DIM>
std::set<unsigned> AbstractTissue<DIM>::GetCellAncestors()
{
    std::set<unsigned> remaining_ancestors;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        remaining_ancestors.insert(cell_iter->GetAncestor());
    }
    return remaining_ancestors;
}

template<unsigned DIM>
c_vector<unsigned, NUM_CELL_MUTATION_STATES> AbstractTissue<DIM>::GetCellMutationStateCount()
{
    return mCellMutationStateCount;
}

template<unsigned DIM>
c_vector<unsigned, NUM_CELL_TYPES> AbstractTissue<DIM>::GetCellTypeCount()
{
    return mCellTypeCount;
}

template<unsigned DIM>
c_vector<unsigned, 5> AbstractTissue<DIM>::GetCellCyclePhaseCount()
{
    return mCellCyclePhaseCount;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::IsGhostNode(unsigned index)
{
    return false;
}

template<unsigned DIM>
TissueCell& AbstractTissue<DIM>::rGetCellUsingLocationIndex(unsigned index)
{
    return *(mLocationCellMap[index]);
}

template<unsigned DIM>
bool AbstractTissue<DIM>::IsCellAssociatedWithAGhostNode(TissueCell cell)
{
    return false;
}

//////////////////////////////////////////////////////////////////////////////
//                             Iterator class                               //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
TissueCell& AbstractTissue<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
TissueCell* AbstractTissue<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return &(*mCellIter);
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::operator!=(const AbstractTissue<DIM>::Iterator& other)
{
    return mCellIter != other.mCellIter;
}

template<unsigned DIM>
typename AbstractTissue<DIM>::Iterator& AbstractTissue<DIM>::Iterator::operator++()
{
    do
    {
        ++mCellIter;
        if (!IsAtEnd())
        {
            mLocationIndex = mCellIter->GetLocationIndex();
        }
    }
    while (!IsAtEnd() && !IsRealCell());

    return (*this);
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::IsRealCell()
{
    return !(    mrTissue.IsCellAssociatedWithAGhostNode(*mCellIter) 
              || mrTissue.IsCellAssociatedWithADeletedNode(*mCellIter)
              || (*this)->IsDead() );
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrTissue.rGetCells().end();
}

template<unsigned DIM>
AbstractTissue<DIM>::Iterator::Iterator(AbstractTissue& rTissue, std::list<TissueCell>::iterator cellIter)
    : mrTissue(rTissue),
      mCellIter(cellIter)
{
    // Make sure the tissue isn't empty
    assert(mrTissue.rGetCells().size() > 0);
    if (!IsAtEnd())
    {
        mLocationIndex = cellIter->GetLocationIndex();
    }
    // Make sure we start at a real cell
    if (mCellIter == mrTissue.rGetCells().begin() && !IsRealCell())
    {
        ++(*this);
    }
}

template<unsigned DIM>
typename AbstractTissue<DIM>::Iterator AbstractTissue<DIM>::Begin()
{
    return Iterator(*this, this->mCells.begin());
}

template<unsigned DIM>
typename AbstractTissue<DIM>::Iterator AbstractTissue<DIM>::End()
{
    return Iterator(*this, this->mCells.end());
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void AbstractTissue<DIM>::WriteMeshToFile(const std::string &rArchiveDirectory, const std::string &rMeshFileName)
{    
}

template<unsigned DIM>
void AbstractTissue<DIM>::CreateOutputFiles(const std::string &rDirectory,
                                            bool rCleanOutputDirectory,
                                            bool outputCellMutationStates,
                                            bool outputCellTypes,
                                            bool outputCellVariables,
                                            bool outputCellCyclePhases,
                                            bool outputCellAncestors)
{
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes");
    mpVizCellTypesFile = output_file_handler.OpenOutputFile("results.vizcelltypes");
    
    if (outputCellAncestors)
    {
        mpCellAncestorsFile = output_file_handler.OpenOutputFile("results.vizancestors");
    }
    if (outputCellMutationStates)
    {
        mpCellMutationStatesFile = output_file_handler.OpenOutputFile("cellmutationstates.dat");
        *mpCellMutationStatesFile << "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
    if (outputCellTypes)
    {
        mpCellTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    }
    if (outputCellVariables)
    {
        mpCellVariablesFile = output_file_handler.OpenOutputFile("cellvariables.dat");
    }
    if (outputCellCyclePhases)
    {
        mpCellCyclePhasesFile = output_file_handler.OpenOutputFile("cellcyclephases.dat");
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::CloseOutputFiles(bool outputCellMutationStates,
                                           bool outputCellTypes,
                                           bool outputCellVariables,
                                           bool outputCellCyclePhases,
                                           bool outputCellAncestors)
{
    mpVizNodesFile->close();
    mpVizCellTypesFile->close();
    
    if (outputCellMutationStates)
    {
        mpCellMutationStatesFile->close();
    }
    if (outputCellTypes)
    {
        mpCellTypesFile->close();
    }
    if (outputCellVariables)
    {
        mpCellVariablesFile->close();
    }
    if (outputCellCyclePhases)
    {
        mpCellCyclePhasesFile->close();
    }
    if (outputCellAncestors)
    {
        mpCellAncestorsFile->close();
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::GenerateCellResults(unsigned locationIndex,
                                              bool outputCellMutationStates,
                                              bool outputCellTypes,
                                              bool outputCellVariables,
                                              bool outputCellCyclePhases,
                                              bool outputCellAncestors,
                                              std::vector<unsigned>& rCellTypeCounter,
                                              std::vector<unsigned>& rCellMutationStateCounter,
                                              std::vector<unsigned>& rCellCyclePhaseCounter)
{
    unsigned colour = STEM_COLOUR;
    if (IsGhostNode(locationIndex) == true)
    {
        colour = INVISIBLE_COLOUR;
    }
    else
    {            
        TissueCell* p_cell = mLocationCellMap[locationIndex];
        
        // Cell cycle phase
        if (outputCellCyclePhases)
        {
            switch (p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase())
            {
                case G_ZERO_PHASE:
                    rCellCyclePhaseCounter[0]++;
                    break;
                case G_ONE_PHASE:
                    rCellCyclePhaseCounter[1]++;
                    break;
                case S_PHASE:
                    rCellCyclePhaseCounter[2]++;
                    break;
                case G_TWO_PHASE:
                    rCellCyclePhaseCounter[3]++;
                    break;
                 case M_PHASE:
                    rCellCyclePhaseCounter[4]++;
                    break;
                default:
                    NEVER_REACHED;
            }
        }
        
        // Cell ancestors
        if (outputCellAncestors)
        {
            colour = p_cell->GetAncestor();
            if (colour == UNSIGNED_UNSET)
            {
                // Set the file to -1 to mark this case.
                colour = 1;
                *mpCellAncestorsFile << "-";
            }
            *mpCellAncestorsFile << colour << " ";
        }
    
        // Set colour dependent on cell type
        switch (p_cell->GetCellType())
        {
            case STEM:
                colour = STEM_COLOUR;
                if (outputCellTypes)
                {
                    rCellTypeCounter[0]++;
                }
                break;
            case TRANSIT:
                colour = TRANSIT_COLOUR;
                if (outputCellTypes)
                {
                    rCellTypeCounter[1]++;
                }
                break;
            case DIFFERENTIATED:
                colour = DIFFERENTIATED_COLOUR;
                if (outputCellTypes)
                {
                    rCellTypeCounter[2]++;
                }
                break;
            case APOPTOTIC:
                colour = APOPTOSIS_COLOUR;
                if (outputCellTypes)
                {
                    rCellTypeCounter[3]++;
                }
                break;
            default:
                NEVER_REACHED;
        }
    
        // Override colours for mutant or labelled cells        
        CellMutationState mutation = p_cell->GetMutationState();
        switch (mutation)
        {
            case HEALTHY:
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[0]++;
                }
                break;
            case APC_ONE_HIT:
                colour = EARLY_CANCER_COLOUR;
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[2]++;
                }
                break;
            case APC_TWO_HIT:
                colour = LATE_CANCER_COLOUR;
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[3]++;
                }
                break;
            case BETA_CATENIN_ONE_HIT:
                colour = LATE_CANCER_COLOUR;
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[4]++;
                }
                break;
            case LABELLED:
                colour = LABELLED_COLOUR;
                if (outputCellMutationStates)
                {
                    rCellMutationStateCounter[1]++;
                }
                break;
            default:
                NEVER_REACHED;
        }
    
        if (p_cell->HasApoptosisBegun())
        {
            // For any type of cell set the colour to this if it is undergoing apoptosis.
            colour = APOPTOSIS_COLOUR;
        }
    
        // Write cell variable data to file if required
        if (outputCellVariables)
        {
            std::vector<double> proteins = p_cell->GetCellCycleModel()->GetProteinConcentrations();

            // Loop over cell positions
            /// \todo Note that the output format of mpCellVariablesFile has changed, it now
            //        outputs the location index (either node index or vertex element index) 
            //        associated with the cell. This is so that this code works for vertex- 
            //        as well as cell-centre-based tissues (see #827).
            for (unsigned i=0; i<DIM; i++)
            {
                *mpCellVariablesFile << locationIndex << " ";
            }
            // Loop over cell variables
            for (unsigned i=0; i<proteins.size(); i++)
            {
                *mpCellVariablesFile << proteins[i] << " ";
            }
        }
    }
    *mpVizCellTypesFile << colour << " ";
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteCellResultsToFiles(bool outputCellMutationStates,
                                                  bool outputCellTypes,
                                                  bool outputCellVariables,
                                                  bool outputCellCyclePhases,
                                                  bool outputCellAncestors,
                                                  std::vector<unsigned>& rCellTypeCounter,
                                                  std::vector<unsigned>& rCellMutationStateCounter,
                                                  std::vector<unsigned>& rCellCyclePhaseCounter)
{
    *mpVizCellTypesFile << "\n";
    
    if (outputCellAncestors)
    {
        *mpCellAncestorsFile << "\n";
    }

    // Write cell mutation state data to file if required
    if (outputCellMutationStates)
    {
        for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
        {
            mCellMutationStateCount[i] = rCellMutationStateCounter[i];
            *mpCellMutationStatesFile << rCellMutationStateCounter[i] << "\t";
        }
        *mpCellMutationStatesFile << "\n";
    }

    // Write cell type data to file if required
    if (outputCellTypes)
    {
        for (unsigned i=0; i<NUM_CELL_TYPES; i++)
        {
            mCellTypeCount[i] = rCellTypeCounter[i];
            *mpCellTypesFile << rCellTypeCounter[i] << "\t";
        }
        *mpCellTypesFile << "\n";
    }

    if (outputCellVariables)
    {
        *mpCellVariablesFile << "\n";
    }

    // Write cell cycle phase data to file if required
    if (outputCellCyclePhases)
    {
        for (unsigned i=0; i<5; i++)
        {
            mCellCyclePhaseCount[i] = rCellCyclePhaseCounter[i];
            *mpCellCyclePhasesFile << rCellCyclePhaseCounter[i] << "\t";
        }
        *mpCellCyclePhasesFile << "\n";
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::WriteTimeAndNodeResultsToFiles(bool outputCellMutationStates,
                                                         bool outputCellTypes,
                                                         bool outputCellVariables,
                                                         bool outputCellCyclePhases,
                                                         bool outputCellAncestors,
                                                         std::vector<unsigned>& rCellTypeCounter,
                                                         std::vector<unsigned>& rCellMutationStateCounter,
                                                         std::vector<unsigned>& rCellCyclePhaseCounter)
{
    // Write current simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetTime();
    
    *mpVizNodesFile << time << "\t";
    *mpVizCellTypesFile << time << "\t";
    
    if (outputCellAncestors)
    {
        *mpCellAncestorsFile << time << "\t";
    }
    if (outputCellMutationStates)
    {
        *mpCellMutationStatesFile << time << "\t";
    }
    if (outputCellTypes)
    {
        *mpCellTypesFile << time << "\t";
    }
    if (outputCellVariables)
    {
        *mpCellVariablesFile << time << "\t";
    }
    if (outputCellCyclePhases)
    {
        *mpCellCyclePhasesFile << time << "\t";
    }
    
    // Set up cell type counter
    rCellTypeCounter.reserve(mCellTypeCount.size());
    for (unsigned i=0; i<NUM_CELL_TYPES; i++)
    {
        rCellTypeCounter[i] = 0;
    }

    // Set up cell mutation state counter
    rCellMutationStateCounter.reserve(mCellMutationStateCount.size());
    for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
    {
        rCellMutationStateCounter[i] = 0;
    }

    // Set up cell cycle phase counter
    rCellCyclePhaseCounter.reserve(5);
    for (unsigned i=0; i<5; i++)
    {
        rCellCyclePhaseCounter[i] = 0;
    }
   
    // Write node data to file
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        // Write node data to file
        if ( !(GetNode(node_index)->IsDeleted()) )
        {
            const c_vector<double,DIM>& position = GetNode(node_index)->rGetLocation();

            for (unsigned i=0; i<DIM; i++)
            {
                *mpVizNodesFile << position[i] << " ";
            }
        }
    }
    *mpVizNodesFile << "\n";
}

namespace boost {
namespace serialization {
template<unsigned DIM>
struct is_abstract<AbstractTissue<DIM> > {
    typedef mpl::bool_<true> type;
    BOOST_STATIC_CONSTANT(bool, value=true);
};
}}

#endif /*ABSTRACTTISSUE_HPP_*/
