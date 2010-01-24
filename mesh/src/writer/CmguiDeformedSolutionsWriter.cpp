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


#include "CmguiDeformedSolutionsWriter.hpp"

template<unsigned DIM>
CmguiDeformedSolutionsWriter<DIM>::CmguiDeformedSolutionsWriter(std::string outputDirectory,
                                                                std::string baseName,
                                                                QuadraticMesh<DIM>& rQuadraticMesh) 
    : CmguiMeshWriter<DIM, DIM>(outputDirectory, baseName),
      mpQuadraticMesh(&rQuadraticMesh),
      mFinalCounter(0)
{
}

template<unsigned DIM>
void CmguiDeformedSolutionsWriter<DIM>::WriteInitialMesh()
{
    std::string saved_base_name = this->mBaseName;
    this->mBaseName = this->mBaseName + "_0";
    this->WriteFilesUsingMesh(*mpQuadraticMesh);
    this->mBaseName = saved_base_name;
}

template<unsigned DIM>
void CmguiDeformedSolutionsWriter<DIM>::WriteDeformationPositions(std::vector<c_vector<double,DIM> >& rDeformedPositions,
                                                                  unsigned counter)
{
    mFinalCounter = counter;
    std::stringstream node_file_name_stringstream;
    node_file_name_stringstream <<  this->mBaseName << "_" << counter << ".exnode";

    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name_stringstream.str());

    // Write the node header
    *p_node_file << "Group name: " << this->mGroupName << "\n";
    switch (DIM)
    {
        case 1:
        {
            *p_node_file << CmguiNodeFileHeader1D;
            break;
        }
        case 2:
        {
            *p_node_file << CmguiNodeFileHeader2D;
            break;
        }
        case 3:
        {
            *p_node_file << CmguiNodeFileHeader3D;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    // Write each node's data
    for (unsigned index=0; index<mpQuadraticMesh->GetNumNodes(); index++)
    {
       *p_node_file << "Node:\t" << index+1 << "\t";

        for (unsigned i=0; i<DIM; i++)
        {
            *p_node_file << rDeformedPositions[index](i) << "\t";
        }
        *p_node_file << "\n";
    }
    p_node_file->close();
}

template<unsigned DIM>
void CmguiDeformedSolutionsWriter<DIM>::WriteCmguiScript()
{
    out_stream p_script_file = this->mpOutputFileHandler->OpenOutputFile("LoadSolutions.com");
    *p_script_file << "#\n# Cmgui script automatically generated by Chaste\n#\n"
                   << "for ($i=0; $i<=" << mFinalCounter << "; $i++) { \n"
                   << "  gfx read node " << this->mBaseName << "_$i time $i\n"
                   << "}\n"
                   << "gfx read ele " << this->mBaseName << "_0\n"
                   << "gfx cr win\n\n";        
    p_script_file->close();
}


template class CmguiDeformedSolutionsWriter<2>;
template class CmguiDeformedSolutionsWriter<3>;
