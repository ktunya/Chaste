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


#ifndef MONODOMAINPDE_HPP_
#define MONODOMAINPDE_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <vector>
#include "AbstractCardiacPde.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "Node.hpp"
#include "Element.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>


/**
 * MonodomainPde class.
 *
 * The monodomain equation is of the form:
 * A (C dV/dt + Iionic) +Istim = Div( sigma_i Grad(V) )
 *
 * where A is the surface area to volume ratio         (1/cm),
 *       C is the capacitance                          (uF/cm^2),
 *       sigma_i is the intracellular conductivity     (mS/cm),
 *       I_ionic is the ionic current                  (uA/cm^2),
 *       I_stim is the intracellular stimulus current  (uA/cm^3).
 *
 * Note that default values of A, C and sigma_i are stored in the parent class
 */
template <unsigned ELEM_DIM, unsigned SPACE_DIM = ELEM_DIM>
class MonodomainPde : public virtual AbstractCardiacPde<ELEM_DIM,SPACE_DIM>, public AbstractLinearParabolicPde<ELEM_DIM, SPACE_DIM>
{
private:
    friend class TestMonodomainPde;

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
    }


public:
    /// Constructor
    MonodomainPde(AbstractCardiacCellFactory<ELEM_DIM,SPACE_DIM>* pCellFactory);
    
    // Another constructor (for archiving)
    MonodomainPde(std::vector<AbstractCardiacCell*> & rCellsDistributed);

    //The following are hidden from the coverage test while it is waiting
    //for a re-factor. (Ticket #157)
#define COVERAGE_IGNORE
    /**
     * This should not be called; use
     * ComputeLinearSourceTermAtNode instead
     */
    double ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& );

    /**
     * This should not be called; use
     * ComputeNonlinearSourceTermAtNode instead
     */
    double ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>& , double );
#undef COVERAGE_IGNORE

   /**
     * Compute the diffusion term at a given point.
     * 
     * @param rX The point in space at which the diffusion term is computed.
     * @param pElement the element for which to compute the contribution
     * @return A matrix.
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(
                const ChastePoint<SPACE_DIM>& rX,
                Element<ELEM_DIM,SPACE_DIM>* pElement);


    double ComputeNonlinearSourceTermAtNode(const Node<SPACE_DIM>& node, double );


    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& );
};

// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS2(MonodomainPde, 1, 1) 
EXPORT_TEMPLATE_CLASS2(MonodomainPde, 2, 2) 
EXPORT_TEMPLATE_CLASS2(MonodomainPde, 3, 3) 

namespace boost
{
namespace serialization
{

template<class Archive, unsigned ELEM_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const MonodomainPde<ELEM_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
 
    const std::vector<AbstractCardiacCell*> & r_cells_distributed = t->GetCellsDistributed();
    
    ar << r_cells_distributed;
    
    /// \todo: #98 Archive intra conductivity tensors       
}    
    
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive, unsigned ELEM_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, MonodomainPde<ELEM_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    std::vector<AbstractCardiacCell*> cells_distributed;

    ar >> cells_distributed;

    /// \todo: #98 Retrieve intra conductivity tensors

    ::new(t)MonodomainPde<ELEM_DIM, SPACE_DIM>(cells_distributed /*, pass intracellular conductivity tensors*/);
}
}
} // namespace ...


#endif /*MONODOMAINPDE_HPP_*/
