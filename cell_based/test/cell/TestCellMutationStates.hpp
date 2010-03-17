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
#ifndef TESTCELLMUTATIONSTATES_HPP_
#define TESTCELLMUTATIONSTATES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "LabelledCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"

class TestCellMutationStates : public AbstractCellBasedTestSuite
{
public:

	void TestCellMutationStateMethods() throw(Exception)
	{
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
		p_state->IncrementCellCount();
		TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
		p_state->DecrementCellCount();
		TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
		TS_ASSERT_THROWS_THIS(p_state->DecrementCellCount(),
				"Cannot decrement cell count: no cells have this mutation state.");
		TS_ASSERT_EQUALS(p_state->GetColour(), 0u);

		TS_ASSERT_EQUALS(p_state->IsType<WildTypeCellMutationState>(), true);
		TS_ASSERT_EQUALS(p_state->IsType<ApcOneHitCellMutationState>(), false);
	}

	void TestArchiveCellMutationState() throw(Exception)
	{
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mutation.arch";

		// Archive a mutation state
		{
			ApcOneHitCellMutationState* p_state = new ApcOneHitCellMutationState();
			p_state->IncrementCellCount();

			TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
			TS_ASSERT_EQUALS(p_state->GetColour(), 0u);

			// Create an output archive
			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			// Write the cell to the archive
			output_arch << p_state;
		}

		// Restore mutation state
		{
			// Initialize a mutation state

			ApcOneHitCellMutationState* p_state;

			// Restore the mutation state
			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive input_arch(ifs);

			input_arch >> p_state;

			TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
			TS_ASSERT_EQUALS(p_state->GetColour(), 0u);

			// Tidy up
			delete p_state;
		}
	}

	void TestArchiveBoostSharedPointer() throw(Exception)
	{
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mutation.arch";

		// Archive a boost::shared_ptr to a mutation state
		{
			boost::shared_ptr<WildTypeCellMutationState> p_state(new WildTypeCellMutationState);
			p_state->IncrementCellCount();

			const boost::shared_ptr<WildTypeCellMutationState> p_const = p_state;

			// Create an output archive
			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			// Write the cell to the archive
			output_arch << p_const;
		}

		// Restore boost::shared_ptr to mutation state
		{
			// Initialize a mutation state
			boost::shared_ptr<WildTypeCellMutationState> p_state;

			// Restore the mutation state
			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive input_arch(ifs);

			input_arch >> p_state;

			TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
			TS_ASSERT_EQUALS(p_state->GetColour(), 0u);
		}
	}

	void TestArchiveBoostSharedPointerUsingAbstract() throw(Exception)
	{
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mutation_abs.arch";

		// Archive a boost::shared_ptr to a mutation state
		{
			boost::shared_ptr<AbstractCellMutationState> p_state(new BetaCateninOneHitCellMutationState);
			p_state->IncrementCellCount();

			const boost::shared_ptr<AbstractCellMutationState> p_const = p_state;

			TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
			TS_ASSERT_EQUALS(p_state->GetColour(), 4u);

			// Create an output archive
			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			// Write the cell to the archive
			output_arch << p_const;
		}

		// Restore boost::shared_ptr to mutation state
		{
			// Initialize a mutation state
			boost::shared_ptr<AbstractCellMutationState> p_state;

			// Restore the mutation state
			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive input_arch(ifs);

			input_arch >> p_state;

			TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
			TS_ASSERT_EQUALS(p_state->GetColour(), 4u);
		}
	}
};

#endif /* TESTCELLMUTATIONSTATES_HPP_ */
