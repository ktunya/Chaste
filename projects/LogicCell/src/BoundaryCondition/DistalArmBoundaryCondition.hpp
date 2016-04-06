/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
* Neither the name of the University of Oxford nor the names of its
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef DISTALARMBOUNDARYCONDITION_HPP_
#define DISTALARMBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
* A boundary condition that confines germ cells to the distal arm of the adult C. elegans germ line.
* On reaching the turn, cells are removed.
*/

template<unsigned DIM>
class DistalArmBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{

private:

	/** Needed for serialization. */
	friend class boost::serialization::access;
	/**
	* Serialize the object.
	* @param archive the archive
	* @param version the current version of this class
	*/
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
		archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
	}

	//Maximum distance that cells can be from the distal arm midline
	double TubeRadius;
	//Length of the tube
	double TubeLength;

public:

	/**
	* Constructor.
	*
	* @param pCellPopulation pointer to the cell population
	* @param radius tube radius
	* @param length tube length
	*/
	DistalArmBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, double radius, double length);


	/**
	* Overridden ImposeBoundaryCondition() method.
	* Apply the cell population boundary condition.
	*
	* @param rOldLocations the node locations before any boundary conditions are applied
	*/
	void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);


	/**
	* Overridden VerifyBoundaryCondition() method.
	* Verify the boundary conditions are consistent.
	* This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
	*
	* @return whether the boundary conditions are satisfied.
	*/
	bool VerifyBoundaryCondition();


	//Getters for private members
	double GetTubeRadius() const;
	double GetTubeLength() const;


	/**
	* Overridden OutputCellPopulationBoundaryConditionParameters() method.
	* Output cell population boundary condition parameters to file.
	*
	* @param rParamsFile the file stream to which the parameters are output
	*/
	void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DistalArmBoundaryCondition)


namespace boost
{
	namespace serialization
	{
		/**
		* Serialize information required to construct a LeaderCellBoundaryCondition.
		*/
		template<class Archive, unsigned DIM>
		inline void save_construct_data(
			Archive & ar, const DistalArmBoundaryCondition<DIM>* t, const BOOST_PFTO unsigned int file_version)
		{
			// Save data required to construct instance
			const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
			ar << p_cell_population;
			const double radius = t->GetTubeRadius();
			ar << radius;
			const double length = t->GetTubeLength();
			ar << length;
		}


		/**
		* De-serialize constructor parameters and initialize a LeaderCellBoundaryCondition.
		*/
		template<class Archive, unsigned DIM>
		inline void load_construct_data(
			Archive & ar, DistalArmBoundaryCondition<DIM>* t, const unsigned int file_version)
		{
			// Retrieve data from archive required to construct new instance
			AbstractCellPopulation<DIM>* p_cell_population;
			ar >> p_cell_population;
			double radius;
			ar >> radius;
			double length;
			ar >> length;

			// Invoke inplace constructor to initialise instance
			::new(t)DistalArmBoundaryCondition<DIM>(p_cell_population, radius, length);
		}
	}
} // namespace ...

#endif /*DISTALARMBOUNDARYCONDITION_HPP_*/