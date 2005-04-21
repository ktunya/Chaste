#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "../TrianglesMeshReader.hpp"

static		AbstractMeshReader *spMeshReader;
class TestTrianglesMeshReaders : public CxxTest::TestSuite
{
	public:
	void testFilesOpen(void)
	{
		
		
		
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements"));
		
	
	}
	
	void testNodesDataRead(void)
	{
		
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumNodes() == 543); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_nodes_disk_522__elements_indexed_from_1"));		
		
	}
	
	void testElementsDataRead(void)
	{
		
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumElements() == 984); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_elements_disk_522_elements_indexed_from_1"));
	
	
	}
	
	void testFacesDataRead(void)
	{
		
		spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements_indexed_from_1");
		
		TS_ASSERT( spMeshReader->GetNumFaces() == 1526); 
		
		
		TS_ASSERT_THROWS_ANYTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/bad_faces_disk_522__elements_indexed_from_1"));		
		
	}
	
	void test3dDataRead(void)
	{
		
			
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/slab_138_elements"));
		
		
		TS_ASSERT (spMeshReader->GetNumElements() == 138);
		
	}
	
	void testIndexFromZero(void)
	{
		
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements"));
		//spMeshReader=new TrianglesMeshReader("pdes/tests/meshdata/disk_522__elements");
		
		TS_ASSERT(spMeshReader->GetMaxNodeIndex() == spMeshReader->GetNumNodes() - 1);
		TS_ASSERT(spMeshReader->GetMinNodeIndex() == 0);
		
	}
	
	
	void testIndexFromOne(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spMeshReader=new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_522_elements_indexed_from_1"));
		
		//spMeshReader=new TrianglesMeshReader("pdes/tests/meshdata/disk_522__elements_indexed_from_1");
		
		TS_ASSERT(spMeshReader->GetMaxNodeIndex() == spMeshReader->GetNumNodes() - 1);
		TS_ASSERT(spMeshReader->GetMinNodeIndex() == 0);
		
	}
	
	
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
