// TestMemfemMeshReader.hpp

/**
 * Test suite for the MemfemMeshReader class.
 * 
 */

#ifndef _TESTMEMFEMMESHREADER_HPP_
#define _TESTMEMFEMMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.hpp"

class TestMemfemMeshReaders : public CxxTest::TestSuite
{
	public:
	
	/**
	 * Check that input files are opened correctly.
	 * 
	 */
	void TestFilesOpen(void)
	{
		MemfemMeshReader *pMeshReader;
		TS_ASSERT_THROWS_NOTHING(
		                  pMeshReader = new MemfemMeshReader(
		                  "pdes/tests/meshdata/Memfem_slab"));		                  
		             
		TS_ASSERT(pMeshReader->GetNumNodes() == 381);
		TS_ASSERT(pMeshReader->GetNumElements() == 1030);
		TS_ASSERT(pMeshReader->GetNumBoundaryFaces() == 758);
		
		std::vector<int> NextBoundaryFace;
		                  
		NextBoundaryFace = pMeshReader->GetNextBoundaryFace();
		
		TS_ASSERT( NextBoundaryFace[0] == 338  );
		TS_ASSERT( NextBoundaryFace[1] == 23 );
		TS_ASSERT( NextBoundaryFace[2] == 374 );
		
		TS_ASSERT(pMeshReader->GetMaxNodeIndex() == pMeshReader->GetNumNodes() - 1);
		
		TS_ASSERT(pMeshReader->GetMinNodeIndex() == 0);			    		
		
		delete pMeshReader;
	}

};

#endif //_TESTMEMFEMMESHREADER_HPP_
