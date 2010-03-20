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
#include "Exception.hpp"
#include "CmguiMeshWriter.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::CmguiMeshWriter(const std::string &rDirectory,
                                                        const std::string &rBaseName,
                                                        const bool &rCleanDirectory)
        : AbstractTetrahedralMeshWriter<ELEMENT_DIM,SPACE_DIM>(rDirectory, rBaseName, rCleanDirectory)
{
    this->mIndexFromZero = false;
    mGroupName = this->mBaseName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::WriteFiles()
{
    //////////////////////////
    // Write the exnode file
    //////////////////////////
    std::string node_file_name = this->mBaseName + ".exnode";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);

    WriteNodeFileHeader(p_node_file);

    // Write each node's data
    for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
    {
        std::vector<double> current_item = this->GetNextNode();

        *p_node_file << "Node:\t" << item_num+1 << "\t";
        for (unsigned i=0; i<ELEMENT_DIM; i++)
        {
            *p_node_file << current_item[i] << "\t";
        }

        *p_node_file << "\n";
    }
    p_node_file->close();

    //////////////////////////
    // Write the exlem file
    //////////////////////////
    std::string elem_file_name = this->mBaseName + ".exelem";
    out_stream p_elem_file = this->mpOutputFileHandler->OpenOutputFile(elem_file_name);

    // Write the elem header
    *p_elem_file << "Group name: " << mGroupName << "\n";
    switch (ELEMENT_DIM)
    {
        case 1:
        {
            *p_elem_file << CmguiElementFileHeader1D;
            break;
        }
        case 2:
        {
            *p_elem_file << CmguiElementFileHeader2D;
            break;
        }
        case 3:
        {
            *p_elem_file << CmguiElementFileHeader3D;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }


    //now we need to figure out how many additional fields we have
    unsigned number_of_fields = mAdditionalFieldNames.size();
    std::stringstream string_of_number_of_fields;
    //we write the number of additional fields + 1 because the coordinates field gets written anyway
    string_of_number_of_fields << number_of_fields+1;
    //and write accordingly the total number of fields
    *p_elem_file << " #Fields="<<string_of_number_of_fields.str()<<"\n";

    //first field (the coordinates field is fixed and always there
    switch (ELEMENT_DIM)
    {
        case 1:
        {
            *p_elem_file << CmguiCoordinatesFileHeader1D;
            break;
        }
        case 2:
        {
            *p_elem_file << CmguiCoordinatesFileHeader2D;
            break;
        }
        case 3:
        {
            *p_elem_file << CmguiCoordinatesFileHeader3D;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }


    //now write the specification for each additional field
    for (unsigned i = 0; i <  number_of_fields; i++)
    {
        //unsigned to string
        std::stringstream i_string;
        i_string << i+2;
        *p_elem_file<<i_string.str()<<")  "<<mAdditionalFieldNames[i]<<" ,";
        switch (ELEMENT_DIM)
        {
            case 1:
            {
                *p_elem_file << CmguiAdditonalFieldHeader1D;
                break;
            }
            case 2:
            {
                *p_elem_file << CmguiAdditonalFieldHeader2D;
                break;
            }
            case 3:
            {
                *p_elem_file << CmguiAdditonalFieldHeader3D;
                break;
            }
            default:
            {
                NEVER_REACHED;
            }
        }

    }

    // Write each elements's data
    for (unsigned item_num=0; item_num<this->GetNumElements(); item_num++)
    {
        std::vector<unsigned> current_element = this->GetNextElement().NodeIndices;

        *p_elem_file << "Element:\t" << item_num+1 << " 0 0 Nodes:\t";
        for (unsigned i=0; i<(ELEMENT_DIM+1); i++)
        {
            *p_elem_file << current_element[i]+1 << "\t";
        }

        *p_elem_file << "\n";
    }
    p_elem_file->close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::SetAdditionalFieldNames(std::vector<std::string>& rFieldNames)
{
    mAdditionalFieldNames = rFieldNames;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::WriteNodeFileHeader(out_stream& rpNodeFile)
{
    // Write the node header
    *rpNodeFile << "Group name: " << this->mGroupName << "\n";
    switch (SPACE_DIM)
    {
        case 1:
        {
            *rpNodeFile << CmguiNodeFileHeader1D;
            break;
        }
        case 2:
        {
            *rpNodeFile << CmguiNodeFileHeader2D;
            break;
        }
        case 3:
        {
            *rpNodeFile << CmguiNodeFileHeader3D;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class CmguiMeshWriter<1,1>;
template class CmguiMeshWriter<1,2>;
template class CmguiMeshWriter<1,3>;
template class CmguiMeshWriter<2,2>;
template class CmguiMeshWriter<2,3>;
template class CmguiMeshWriter<3,3>;
