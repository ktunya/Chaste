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

#include <algorithm>

#include "CellMutationStateRegistry.hpp"

CellMutationStateRegistry* CellMutationStateRegistry::mpInstance = NULL;

CellMutationStateRegistry* CellMutationStateRegistry::Instance()
{
	if (mpInstance == NULL)
	{
		mpInstance = new CellMutationStateRegistry;
	}
	return mpInstance;
}

const std::vector<boost::shared_ptr<AbstractCellMutationState> >& CellMutationStateRegistry::rGetAllMutationStates()
{
	return mMutationStates;
}

void CellMutationStateRegistry::Clear()
{
	mMutationStates.clear();
	mOrderingHasBeenSpecified = false;
}

CellMutationStateRegistry::CellMutationStateRegistry()
	: mOrderingHasBeenSpecified(false)
{
}


CellMutationStateRegistry* CellMutationStateRegistry::TakeOwnership()
{
	mpInstance = NULL;
	return this;
}

void CellMutationStateRegistry::SpecifyOrdering(const std::vector<boost::shared_ptr<AbstractCellMutationState> >& rOrdering)
{
	if (mOrderingHasBeenSpecified)
	{
		EXCEPTION("An ordering has already been specified.");
	}
	for (unsigned i=0; i<mMutationStates.size(); i++)
	{
		std::vector<boost::shared_ptr<AbstractCellMutationState> >::const_iterator it
			= find(rOrdering.begin(), rOrdering.end(), mMutationStates[i]);
		if (it == rOrdering.end())
		{
			EXCEPTION("The given ordering doesn't include all mutation states in the registry.");
		}
	}
	mMutationStates = rOrdering;

	mOrderingHasBeenSpecified = true;
}

bool CellMutationStateRegistry::HasOrderingBeenSpecified()
{
	return mOrderingHasBeenSpecified;
}

