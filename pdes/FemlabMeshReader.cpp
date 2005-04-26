/** Implementation file for the FemlabMeshReader class.*/

#include "FemlabMeshReader.hpp"
/**
 * The constructor takes the path to and names of a set of Femlab mesh files 
 * (ie. the node, elements and face files (in that order) and allows the data to
 * be queried.
 * Typical use:
 *    AbstractMeshReader *spMeshReader=new FemlabMeshReader(
 * 						  "pdes/tests/meshdata/",
 *		                  "femlab_lshape_nodes.dat",
 * 						  "femlab_lshape_elements.dat",
 * 						  "femlab_lshape_edges.dat",);
 * Also calls the superclass AbstractMeshReader's constructor
 */ 
FemlabMeshReader::FemlabMeshReader(std::string pathBaseName, std::string nodeFileName, std::string elementFileName, std::string edgeFileName)
{
	bool indexed_from_zero = false;
	bool already_checked_indexing = false;
	
	int num_attributes, max_marker;
	
	//Copy path and base name of files to private data
	mPathBaseName=pathBaseName;

	//Open node file and store the lines as a vector of strings (minus the comments) 	
	nodeFileName=pathBaseName+nodeFileName;
	mNodeRawData=GetRawDataFromFile(nodeFileName);
	
	// Read the node data using TokenizeStringsToDoubles method
	mNodeData = TokenizeStringsToDoubles(mNodeRawData);
	mNumNodes = mNodeData.size(); 
	//Initialise iterator for public GetNextNode method
	mpNodeIterator = mNodeData.begin();
	
	/*//Check that the size of the data matches the information in the header
	if (mNumNodes != mNodeData.size())
	{
		throw Exception("Number of nodes does not match expected number declared in header");
	}*/
	
	//Open element file and store the lines as a vector of strings (minus the comments) 	
	elementFileName=pathBaseName+elementFileName;
	mElementRawData=GetRawDataFromFile(elementFileName);

 	/*//Only order-1 triangles or tetrahedra are currently supported
	if (mNumElementNodes != mDimension+1)
	{
		throw Exception("Number of nodes per element is not supported");
	}*/
	
	// Read the rest of the element data using TokenizeStringsToInts method
	mElementData = TokenizeStringsToInts(mElementRawData,mDimension+1);
	mNumElements = mElementData.size();
 	mpElementIterator = mElementData.begin();
 	
 	
 	/*//Check that the size of the data matches the information in the header
 	if (mNumElements != mElementData.size())
	{
		throw Exception("Number of elements does not match expected number declared in header");
	}*/

	/*Open face file and store the lines as a vector of strings (minus the comments) 	
	 * Note that the Triangles equivalent of a face file is an edge file.
	 * The two file formats are similar.  We store edges as "faces" but the superclass
	 * provides a GetNextEdge method which queries this data.
	 */

	edgeFileName = pathBaseName+edgeFileName;
	mFaceRawData=GetRawDataFromFile(edgeFileName);
	
	// Read the rest of the face/edge data using TokenizeStringsToInts method
	mFaceData = TokenizeStringsToInts(mFaceRawData,mDimension);
	mNumFaces = mFaceData.size();	
	mpFaceIterator = mFaceData.begin();
	
	mBoundaryFaceData = mFaceData; // Femlab returns only the boundary faces
	mNumBoundaryFaces = mBoundaryFaceData.size();
	mpBoundaryFaceIterator = mBoundaryFaceData.begin();
	
	/*//Check that the size of the data matches the information in the header
	if (mNumFaces != mFaceData.size())
	{
		throw Exception("Number of faces does not match expected number declared in header");
	}*/
	
	
	
}

/** 
 * TokenizeStringsToDoubles is specific to reading node data which came from 
 * a Triangles or Tetgen file.
 * http://tetgen.berlios.de/fformats.node.html
 * http://www-2.cs.cmu.edu/~quake/triangle.node.html
 *  Each string is expected to be an index number, 2 or 3 doubles (representing x,y,z)
 *  an optional list of attributes and an optional boundary marker
 * 	NB: Attributes and boundary markers are currently ignored.
 * Return value is a vector where each item is a vector of double which represents 
 * position.  Indices are implicit in the vector.
 */
	
std::vector<std::vector<double> > FemlabMeshReader::TokenizeStringsToDoubles(
												std::vector<std::string> rawData)
{
	std::vector<std::vector<double> > tokenized_data; // Output
	
	//Iterate over the lines of input
	int dimension_count = 0;
    std::vector<std::string>::iterator the_iterator;
   	for( the_iterator = rawData.begin(); the_iterator != rawData.end(); the_iterator++ )
   	{
   		std::string line_of_data=*the_iterator;
     	std::stringstream line_stream(line_of_data);
     	
     	if (dimension_count == 0)
   		{
   			//First iteration, build the tokenized_data vector and push in x coordinates
			while (!line_stream.eof())	   		
     		{
     			double item_coord;
     			
				std::vector<double> x_coord; 
     			line_stream >> item_coord;
     			x_coord.push_back(item_coord);
     			tokenized_data.push_back(x_coord);
     		}
   		}
   		else
   		{
     		int current_node = 0;
     		
     		//Other iterations, push in coordinates other than x.
			while (!line_stream.eof())
     		{
				double item_coord; 
     			line_stream >> item_coord;
     			tokenized_data[current_node].push_back(item_coord);
     			current_node++;
     		}
   		}
   		//dimension of mesh is the same as the line of rawData.
     	dimension_count++;
     	     	
    }
    
    mDimension = dimension_count;
    return tokenized_data;
}


/** 
 * TokenizeStringsToInts is for reading element, face or edge data which came from 
 * a Triangles or Tetgen file.
 * http://tetgen.berlios.de/fformats.ele.html
 * http://www-2.cs.cmu.edu/~quake/triangle.ele.html
 * http://tetgen.berlios.de/fformats.face.html
 * http://www-2.cs.cmu.edu/~quake/triangle.edge.html
 *  Each string is expected to be:
 *  an index number, 
 *  2, 3 or 4 node indices
 *  ( In 2-D: 2 indices for an edge, 3 for a triangle)
 *  ( In 3-D: 3 indices for a face, 4 for a tetrahedron)
 *  an optional list of attributes (if it's an element)
 *  and an optional boundary marker (if it's an edge or face)
 * 	NB: Attributes and boundary markers are currently ignored.
 * Return value is a vector where each item is a vector of ints which represents 
 * indices of nodes.  
 */
std::vector<std::vector<int> > FemlabMeshReader::TokenizeStringsToInts(
												std::vector<std::string> rawData,
												int dimensionOfObject)
{
	std::vector<std::vector<int> > tokenized_data;

	/* There are dimensionOfObject lines to be read */	
   	for( int i=0;  i<dimensionOfObject; i++ )
   	{
   		std::string line_of_data=rawData[i];
     	std::stringstream line_stream(line_of_data);
     	
     	if (i == 0)
   		{
   			//Form the vector which represents the position of this item     	     		   		
			while (!line_stream.eof())	   		
     		{
     			double item_index;
     			
				std::vector<int> first_index;
     			line_stream >> item_index;
     			first_index.push_back((int)(item_index-0.5)); //item indices should be minus 1.
     			tokenized_data.push_back(first_index);   			     		
     		}
   		}
   		else
   		{
     		int current_node = 0;
     		
     		//Form the vector which represents the position of this item     	     		   		
			while (!line_stream.eof())
     		{
				double item_index;
     			line_stream >> item_index;
     			tokenized_data[current_node].push_back((int)(item_index-0.5));
     			current_node++;
     		}
   		}   	     	
    }	    
    return tokenized_data;
}


/**
 * Destructor
 */
FemlabMeshReader::~FemlabMeshReader()
{
}
