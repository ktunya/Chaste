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

#include "DistalArmBoundaryCondition.hpp"
#include "NodeBasedCellPopulation.hpp"


//Constructor
template<unsigned DIM>
DistalArmBoundaryCondition<DIM>::DistalArmBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
	double radius,
	double length)
	: AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
	TubeRadius(radius),
	TubeLength(length){

	if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation) == NULL)
	{
		EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition.");
	}
	if (DIM != 3)
	{
		EXCEPTION("This boundary condition is intended for use in 3D.");
	}
}



//Takes cells lying outside the gonad boundary and places them back in the monolayer.
template<unsigned DIM>
void DistalArmBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{

	//Iterate over the cell population
	for (typename AbstractCellPopulation<DIM>::RealCellsIterator cell_iter = this->mpCellPopulation->Begin();
		cell_iter != this->mpCellPopulation->End();
		++cell_iter)
	{

		//Read in some properties of this cell
		Node<DIM>* cell_centre_node = this->mpCellPopulation->GetNode(this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter));
		
		c_vector<double, DIM> cell_location = cell_centre_node->rGetLocation();
		
		double cellradius = cell_centre_node->GetRadius();

		//Identify the closest point on leader cell path to this cell
		if (cell_location[0] > TubeLength){
			//Do nothing. This cell will be marked for death by the killer.	
		}
		else if (cell_location[0] < 0){
			//In the endcap
			double distance = sqrt(cell_location[0] * cell_location[0] + cell_location[1] * cell_location[1] + cell_location[2] * cell_location[2]);
			cell_centre_node->rGetModifiableLocation() = cell_location * ((TubeRadius - cellradius) / distance);
		}
		else{
			//In the middle of the tube
			double distance = sqrt(cell_location[1] * cell_location[1] + cell_location[2] * cell_location[2]);
			c_vector<double, DIM> newLoc;
			newLoc[0] = cell_location[0];
			newLoc[1] = cell_location[1] * ((TubeRadius - cellradius) / distance);
			newLoc[2] = cell_location[2] * ((TubeRadius - cellradius) / distance);
			cell_centre_node->rGetModifiableLocation() = newLoc;
		}

		cell_iter->GetCellData()->SetItem("DistanceAwayFromDTC", TubeRadius + cell_location[0]);
		cell_iter->GetCellData()->SetItem("MaxRadius", TubeRadius);

	}
}




template<unsigned DIM>
double DistalArmBoundaryCondition<DIM>::GetTubeRadius() const{
	return TubeRadius;
};
template<unsigned DIM>
double DistalArmBoundaryCondition<DIM>::GetTubeLength() const{
	return TubeLength;
};



//Boundary condition verification. Required when a BC is used in combination with other, possibly 
//incompatible BCs.
template<unsigned DIM>
bool DistalArmBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
	bool condition_satisfied = true;
	//double maxDistOutside = 0.001;
	//
	//for (typename AbstractCellPopulation<DIM>::RealCellsIterator cell_iter = this->mpCellPopulation->Begin();
	//	cell_iter != this->mpCellPopulation->End();
	//	++cell_iter)
	//{
	//
	//	//Read in some properties of this cell
	//	Node<DIM>* cell_centre_node = this->mpCellPopulation->GetNode(this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter));
	//	c_vector<double, DIM> cell_location = cell_centre_node->rGetLocation();
	//	double cellradius = cell_centre_node->GetRadius();
	//
	//	//Identify the closest point on leader cell path to this cell
	//	if (cell_location[0] > TubeLength){
	//	}
	//	else if (cell_location[0] < 0){
	//		double distance = sqrt(cell_location[0] * cell_location[0] + cell_location[1] * cell_location[1] + cell_location[2] * cell_location[2]);
	//		if (distance > (TubeRadius - cellradius) + maxDistOutside){
	//			std::cout << "Endcap Cell" << "\t" << cell_location[0] << "\t" << cell_location[1] << "\t" << cell_location[2] << std::endl;
	//			condition_satisfied = false;
	//			break;
	//		}
	//	}
	//	else{
	//		double distance = sqrt(cell_location[1] * cell_location[1] + cell_location[2] * cell_location[2]);
	//		if (distance > (TubeRadius - cellradius) + maxDistOutside){
	//			std::cout << "Tube Cell" << "\t" << cell_location[0] << "\t" << cell_location[1] << "\t" << cell_location[2] << std::endl;
	//			condition_satisfied = false;
	//			break;
	//		}
	//	}
	//}

	return condition_satisfied;
}



//Parameter output to log file
template<unsigned DIM>
void DistalArmBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{

	*rParamsFile << "\t\t\t<EndOfTube1>" << 0 << "," << 0 << "," << 0 << "</EndOfTube1>\n";
	*rParamsFile << "\t\t\t<EndOfTube2>" << TubeLength << "," << 0 << "," << 0 << "</EndOfTube2>\n";
	*rParamsFile << "\t\t\t<RadiusOfTube>" << TubeRadius << "</RadiusOfTube>\n";

	// Call method on parent class
	AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DistalArmBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DistalArmBoundaryCondition)