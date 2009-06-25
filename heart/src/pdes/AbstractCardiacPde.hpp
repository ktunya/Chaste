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
#ifndef ABSTRACTCARDIACPDE_HPP_
#define ABSTRACTCARDIACPDE_HPP_

#include <vector>
#include <climits> // Work around a boost bug - see #1024.
#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractConductivityTensors.hpp"

#include "ReplicatableVector.hpp"
#include "HeartConfig.hpp"
#include "ArchiveLocationInfo.hpp"

// Needs to be included last
#include "TemplatedExport.hpp"


//// OLD NOTE: read this if AbstractPde is brought back
// IMPORTANT NOTE: the inheritance of AbstractPde has to be 'virtual'
// ie "class AbstractCardiacPde : public virtual AbstractPde"
// because AbstractPde will be the top class in a 'dreaded diamond':
//      A
//     / \     A = AbstractPde, B = AbstractCardiac, C = AbstractLinearParabolic (etc)
//    B   C    D = MonodomainPde
//     \ /
//      D
//
// B and C must use virtual inheritence of A in order for D to only contain 1 instance
// of the member variables in A

/**
 * Class containing common functionality to monodomain and bidomain PDEs.
 *
 * \todo more details...
 */
template <unsigned ELEM_DIM, unsigned SPACE_DIM = ELEM_DIM>
class AbstractCardiacPde
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //Archive vector factory?  Do we have one available?  (If not - do it in the test to guard against errors)

        //Archive cells using handler.GetProcessUniqueFilePath("") and ArchiveLocationInfo

        archive & mpDistributedVectorFactory;
    }


protected:

    /** Intracellular conductivity tensors.  Not archived, since it's loaded from the HeartConfig singleton. */
    AbstractConductivityTensors<SPACE_DIM> *mpIntracellularConductivityTensors;

    /** The vector of cells. Distributed. */
    std::vector< AbstractCardiacCell* > mCellsDistributed;

    /**
     *  Cache containing all the ionic currents for each node,
     *  replicated over all processes.
     */
    ReplicatableVector mIionicCacheReplicated;
    /**
     *  Cache containing all the stimulus currents for each node,
     *  replicated over all processes.
     */
    ReplicatableVector mIntracellularStimulusCacheReplicated;

    /**
     *  Constant set to 1 in monodomain and 2 in bidomain. Used when accessing
     *  the voltage components in the solution vector (because the solution vector
     *  is of the form (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N), where V_j is
     *  the voltage at node j and phi_j is the extracellular potential at node j.
     */
    const unsigned mStride;

    /** Local pointer to the HeartConfig singleton instance, for convenience. */
    HeartConfig* mpConfig;

    /**
     * Whether we need to replicate the caches.
     *
     * When doing matrix-based RHS assembly, we only actually need information from
     * cells/nodes local to the processor, so replicating the caches is an
     * unnecessary communication overhead.
     *
     * Defaults to true.
     */
    bool mDoCacheReplication;
    /**
     * This is to mark the conventional assembly on the first time step.
     *
     * \todo maybe we don't want the conventional assembly even in the first time step.
     */
    bool mDoOneCacheReplication;

    /**
     * Local pointer to the distributed vector factory associated with the mesh object used.
     *
     * Used to retrieve node ownership range when needed.
     */
    DistributedVectorFactory* mpDistributedVectorFactory;
public:
    /**
     * This constructor is called from the Initialise() method of the CardiacProblem class.
     * It creates all the cell objects, and sets up the conductivities.
     *
     * \todo tidy up using extract method refactoring
     *
     * @param pCellFactory  factory to use to create cells.
     * @param stride  determines how to access \f$V_m\f$ in the solution vector (1 for monodomain, 2 for bidomain).
     */
    AbstractCardiacPde(AbstractCardiacCellFactory<ELEM_DIM,SPACE_DIM>* pCellFactory,
                       const unsigned stride=1);

    /**
     * This constructor is called by the archiver
     *
     * @param rCellsDistributed  pointers to the cardiac cells.
     * @param stride  determines how to access \f$V_m\f$ in the solution vector (1 for monodomain, 2 for bidomain).
     */
    AbstractCardiacPde(std::vector<AbstractCardiacCell*> & rCellsDistributed,
                       const unsigned stride = 1);

    /** Virtual destructor */
    virtual ~AbstractCardiacPde();

    /**
     * Set whether or not to replicate the caches across all processors.
     *
     * See also mDoCacheReplication.
     * @param doCacheReplication - true if the cache needs to be replicated
     */
    void SetCacheReplication(bool doCacheReplication);

    /** Get the intracellular conductivity tensor for the given element
     * @param elementIndex  index of the element of interest
     */
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetIntracellularConductivityTensor(unsigned elementIndex);

    /**
     * Get a pointer to a cell, indexed by the global node index.
     *
     * \note Should only called by the process owning the cell -
     * triggers an assertion otherwise.
     *
     * @param globalIndex  global node index for which to retrieve a cell
     */
    AbstractCardiacCell* GetCardiacCell( unsigned globalIndex );


    /**
     *  SolveCellSystems()
     *
     *  Integrate the cell ODEs and update ionic current etc for each of the
     *  cells, between the two times provided.
     *
     *  \note This used to be called PrepareForAssembleSystem, but
     *  that method is now a virtual method in the assemblers not the
     *  pdes.
     *
     * @param currentSolution  the current voltage solution vector
     * @param currentTime  the current simulation time
     * @param nextTime  when to simulate the cells until
     */
    virtual void SolveCellSystems(Vec currentSolution, double currentTime, double nextTime);

    /** Get the entire ionic current cache */
    ReplicatableVector& rGetIionicCacheReplicated();

    /** Get the entire stimulus current cache */
    ReplicatableVector& rGetIntracellularStimulusCacheReplicated();


    /**
     * Update the Iionic and intracellular stimulus caches.
     *
     * @param globalIndex  global index of the entry to update
     * @param localIndex  local index of the entry to update
     * @param nextTime  the next PDE time point, at which to evaluate the stimulus current
     */
    void UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime);

    /**
     *  Replicate the Iionic and intracellular stimulus caches.
     */
    void ReplicateCaches();

    /**
     *  Returns a reference to the vector of distributed cells. Needed for archiving.
     */
    const std::vector<AbstractCardiacCell*>& GetCellsDistributed() const;

};

// Declare identifier for the serializer
EXPORT_ABSTRACT_TEMPLATE_CLASS_SAME_DIMS(AbstractCardiacPde);
EXPORT_ABSTRACT_TEMPLATE_CLASS2(AbstractCardiacPde, 1, 2);
EXPORT_ABSTRACT_TEMPLATE_CLASS2(AbstractCardiacPde, 1, 3);

#endif /*ABSTRACTCARDIACPDE_HPP_*/

