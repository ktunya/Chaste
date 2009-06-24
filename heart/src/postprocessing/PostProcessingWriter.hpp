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

#ifndef POSTPROCESSINGWRITER_HPP_
#define POSTPROCESSINGWRITER_HPP_

#include "PetscTools.hpp"
#include "Hdf5DataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "HeartConfig.hpp"
#include "Debug.hpp"

#include <iostream>

/** 
 * Write out physiological parameters at the end of a simulation
 * - APD map
 * - ...
 */
class PostProcessingWriter
{
  private:
    PropagationPropertiesCalculator* mpCalculator; /**< PropagationPropertiesCalculator based on HDF5 data reader*/
    unsigned mNumberOfNodes; /**< Number of nodes in the mesh (got from the data reader)*/
    
  public:
    /** 
     * Constructor
     * @param pDataReader  an HDF5 reader from which to build the PropagationPropertiesCalculator
     */
    PostProcessingWriter(Hdf5DataReader* pDataReader)
    {
         mpCalculator = new PropagationPropertiesCalculator(pDataReader);
         mNumberOfNodes = pDataReader->GetNumberOfRows();
    }
    
    /**
     * Destructor
     */
    ~PostProcessingWriter()
    {
        delete mpCalculator;
    }
    
    /**
     * Method for opening an APD map file and writing one row per node
     * line 1: <first APD for node 0> <second APD for node 0> ...
     * line 2: <first APD for node 1> <second APD for node 1> ...
     * etc.
     * 
     * Nodes where there is no APD are respresented by a single
     * 0
     * 
     * @param  threshold - Vm used to signify the upstroke (mV)
     * @param  repolarisationPercentage eg. 90.0 for APD90
     */
    void WriteApdMapFile(double threshold, double repolarisationPercentage)
    {
        if(PetscTools::AmMaster())
        {
            out_stream p_file=out_stream(NULL);
            OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
            std::stringstream stream;
            stream << repolarisationPercentage;
            p_file = output_file_handler.OpenOutputFile("Apd" + stream.str() + "Map.dat");
            for (unsigned node_index = 0; node_index < mNumberOfNodes; node_index++)
            { 
                std::vector<double> apds;
                try
                {
                    apds = mpCalculator->CalculateAllActionPotentialDurations(repolarisationPercentage, node_index, threshold);
                    assert(apds.size() != 0);
                }
                catch(Exception& e)
                {                    
                    apds.push_back(0);
                    assert(apds.size() == 1);
                }
                for (unsigned i = 0; i < apds.size(); i++)
                {
                    *p_file << apds[i] << "\t";
                }
                *p_file << std::endl;
            }
            p_file->close();
        }
    }
    

    /**
     * Write out times of each upstroke for each node:
     * 
     * line 1: <first upstroke time for node 0> <second upstroke time for node 0> ...
     * line 2: <first upstroke time for node 1> <second upstroke time for node 1> ...
     * etc.
     * 
     * If there is no upstroke then there will a ...///\todo Fix (see below)
     * 
     * @param threshold  - Vm used to signify the upstroke (mV) 
     */
    void WriteUpstrokeTimeMap(double threshold)
    {
        if(PetscTools::AmMaster())
        {
            out_stream p_file=out_stream(NULL);
            OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
            p_file = output_file_handler.OpenOutputFile("UpstrokeTimeMap.dat");
            for (unsigned node_index = 0; node_index < mNumberOfNodes; node_index++)
            { 
                std::vector<double> upstroke_times;
                upstroke_times = mpCalculator->CalculateUpstrokeTimes(node_index, threshold);
                assert(upstroke_times.size()!=0); ///\todo Fix (see above)
                for (unsigned i = 0; i < upstroke_times.size(); i++)
                {
                    *p_file << upstroke_times[i] << "\t";
                }
                *p_file << std::endl;
            }
            p_file->close();
        }
    }

    /**
     * Write out velocities of each max upstroke for each node:
     * 
     * line 1: <first upstroke velocity for node 0> <second upstroke velocity for node 0> ...
     * line 2: <first upstroke velocity for node 1> <second upstroke velocity for node 1> ...
     * etc.
     * 
     * If there is no upstroke then there will a ...///\todo Fix (see below)
     * 
     * @param threshold  - Vm used to signify the upstroke (mV) 
     */
    void WriteMaxUpstrokeVelocityMap(double threshold)
    {
        if(PetscTools::AmMaster())
        {
            out_stream p_file=out_stream(NULL);
            OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
            p_file = output_file_handler.OpenOutputFile("MaxUpstrokeVelocityMap.dat");
            for (unsigned node_index = 0; node_index < mNumberOfNodes; node_index++)
            { 
                std::vector<double> upstroke_velocities;
                upstroke_velocities = mpCalculator->CalculateAllMaximumUpstrokeVelocities(node_index, threshold);
                assert(upstroke_velocities.size()!=0); ///\todo Fix (see above)
                for (unsigned i = 0; i < upstroke_velocities.size(); i++)
                {
                    *p_file << upstroke_velocities[i] << "\t";
                }
                *p_file << std::endl;
            }
            p_file->close();
        }
    }
};

#endif /*POSTPROCESSINGWRITER_HPP_*/
