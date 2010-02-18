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
#ifndef TESTHEARTGEOMETRYINFORMATION_HPP_
#define TESTHEARTGEOMETRYINFORMATION_HPP_

#include "TrianglesMeshReader.hpp"
#include "HeartGeometryInformation.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"  //temporary
//#include "MeshalyzerMeshWriter.hpp" //temporary
#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp" //temporary

class TestHeartGeometryInformation : public CxxTest::TestSuite
{
private:
    /* Helper method to write out 'fake' face files given the vectors containing the nodes at the surface*/
    void WriteFakeFaceFile(const std::string& rOutputDir, const std::string& rFilename,
                           const std::vector<unsigned>& rNodeLayers, unsigned spaceDim)
    {
        OutputFileHandler handler(rOutputDir, false);
        
        if (PetscTools::AmMaster())
        {
            out_stream p_file = handler.OpenOutputFile(rFilename);
            p_file->close();
        }
        //Each process may have a small number of the nodes
        for (unsigned proc_turn=0; proc_turn<PetscTools::GetNumProcs(); proc_turn++)
        {    
            if (PetscTools::GetMyRank()==proc_turn)
            {
                out_stream p_file = handler.OpenOutputFile(rFilename, std::ios::app);
                
                for (unsigned i=0; i<rNodeLayers.size(); i++)
                {
                    for (unsigned j=0; j<spaceDim; j++)
                    {
                       * p_file << (j == 0 ? "" : "  ") << (1+rNodeLayers[i]);
                    }
                   * p_file << std::endl;
                }
                p_file->close();
            }
            PetscTools::Barrier();
        }
    }

public:
    void TestCalculateRelativeWallPositionSimple2dMesh() throw(Exception)
    {   
        ParallelTetrahedralMesh<2,2> mesh;
        //This mesh will have 6 nodes per face, spaced by 1
        mesh.ConstructRectangularMesh(5, 5);

        std::vector<unsigned> left_face;
        std::vector<unsigned> right_face;
        
        unsigned low_index=mesh.GetDistributedVectorFactory()->GetLow();
        unsigned high_index=mesh.GetDistributedVectorFactory()->GetHigh();
        
        for (unsigned index=low_index; index<high_index; index++)
        {  
            // Get the nodes at the left face of the square
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
            {
                left_face.push_back(index);
            }
            // Get the nodes at the right face of the square
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                right_face.push_back(index);
            }
            
        }           
        // Write our fake face files
        std::string output_dir = "HeartGeom2d";
        WriteFakeFaceFile(output_dir, "epi.tri", left_face, 2u);
        WriteFakeFaceFile(output_dir, "endo.tri", right_face, 2u);

        PetscTools::Barrier(); // Make sure files are written

        // Read in
        OutputFileHandler handler(output_dir, false);
        std::string dir_path = handler.GetOutputDirectoryFullPath();
        //call the constructor that takes in the surface files...
        HeartGeometryInformation<2> info(mesh, dir_path + "/epi.tri", dir_path + "/endo.tri", false);

        
        //and then we test the method to evaluate the position in the wall (again)
        for (unsigned index=low_index; index<high_index; index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(5.0-x)/5.0);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEpicardium()[index],x);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEndocardium()[index],(5.0-x));
        } 
    }

    void TestCalculateRelativeWallPositionSimple3dMesh() throw(Exception)
    {
        TetrahedralMesh<3,3> mesh;
        //This mesh will have 6 nodes per face, spaced by 1
        mesh.ConstructCuboid(5, 5, 5);

        std::vector<unsigned> left_face;
        std::vector<unsigned> right_face;

        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {  
            // Get the nodes at the left face of the cube
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
            {
                left_face.push_back(index);
            }
            // Get the nodes at the right face of the cube
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                right_face.push_back(index);
            }
            
        }
        std::string output_dir = "HeartGeom3d";
        WriteFakeFaceFile(output_dir, "epi.tri", left_face, 3u);
        WriteFakeFaceFile(output_dir, "endo.tri", right_face, 3u);          
        OutputFileHandler handler(output_dir, false);
        std::string dir_path = handler.GetOutputDirectoryFullPath();
        //call the constructor that takes in the surface files...
        HeartGeometryInformation<3> info(mesh, dir_path + "/epi.tri", dir_path + "/endo.tri", false);
        
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(5-x)/5);
        } 
    }
    
    void TestCalculateRelativeWallPositionWithSurfaces() throw(Exception)
    {
        TetrahedralMesh<3,3> mesh;
        //This mesh will have 9 nodes per side, spaced by 1, it is a cube
        mesh.ConstructCuboid(8, 8, 8);

        std::vector<unsigned> epi_face;
        std::vector<unsigned> lv_face;
        std::vector<unsigned> rv_face;
        //Define three surfaces, epi, lv and rv.
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {  
            // Get the nodes at cube face considered to be epi (at both external faces)
            if (  (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
                ||(fabs(mesh.GetNode(index)->rGetLocation()[0]-8.0) < 1e-6))
            {
                epi_face.push_back(index);
            }
            // Get the nodes at cube face considered to be lv (at the plane defined by x=3)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-3.0) < 1e-6)
            {
                lv_face.push_back(index);
            }
            // Get the nodes at cube face considered to be rv (at the plane defined by x=5)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                rv_face.push_back(index);
            }
            
        }           
        
        
        // Write our fake face files
        std::string output_dir = "HeartGeom3d";
        WriteFakeFaceFile(output_dir, "epi.tri", epi_face, 3u);
        WriteFakeFaceFile(output_dir, "lv.tri", lv_face, 3u);
        WriteFakeFaceFile(output_dir, "rv.tri", rv_face, 3u);

        PetscTools::Barrier(); // Make sure files are written

        // Read in
        OutputFileHandler handler(output_dir, false);
        std::string dir_path = handler.GetOutputDirectoryFullPath();
        HeartGeometryInformation<3> info(mesh, dir_path + "/epi.tri", dir_path + "/lv.tri", dir_path + "/rv.tri", false);
        
        //first we test the get methods for the nodes on the surfaces
        std::vector<unsigned> nodes_on_lv = info.rGetNodesOnLVSurface();
        std::vector<unsigned> nodes_on_rv = info.rGetNodesOnRVSurface();
        std::vector<unsigned> nodes_on_epi = info.rGetNodesOnEpiSurface();
        TS_ASSERT_EQUALS(nodes_on_lv.size(),lv_face.size());
        TS_ASSERT_EQUALS(nodes_on_rv.size(),rv_face.size());
        TS_ASSERT_EQUALS(nodes_on_epi.size(),epi_face.size());
        //check that the vectors filled in by the constructor are the same as the original ones
        for (unsigned i = 0; i < lv_face.size();i++)
        {
            TS_ASSERT_EQUALS(nodes_on_lv[i],lv_face[i]);
            TS_ASSERT_EQUALS(nodes_on_rv[i],rv_face[i]);
            TS_ASSERT_EQUALS(nodes_on_epi[i],epi_face[i]);
        }

        //and then we test the method to evaluate the position in the wall (again)
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            //in the lv...
            if (x<=3)
            {
                TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(3.0-x)/3.0);
            }
            //..in the septum...
            else if ((x>3)&&(x<5))
            {
                TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index), 1.0/5.0);
            }
            //...and in the rv.
            else if (x>=5)
            {
                TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(x-5.0)/3.0);
            }
        } 
    }
    
    void TestDetermineLayerForEachNodeWritingAndReading() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        //This mesh will have 31 nodes per side, spaced by 1, it is a cube
        mesh.ConstructCuboid(30, 30, 30);

        std::vector<unsigned> epi_face;
        std::vector<unsigned> lv_face;
        std::vector<unsigned> rv_face;
        //Define three surfaces, epi, lv and rv.
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {  
            // Get the nodes at cube face considered to be epi (at both external faces)
            if (  (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
                ||(fabs(mesh.GetNode(index)->rGetLocation()[0]-30.0) < 1e-6))
            {
                epi_face.push_back(index);
            }
            // Get the nodes at cube face considered to be lv (at the plane defined by x=10)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-10.0) < 1e-6)
            {
                lv_face.push_back(index);
            }
            // Get the nodes at cube face considered to be rv (at the plane defined by x=20)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-20.0) < 1e-6)
            {
                rv_face.push_back(index);
            }
            
        }  
        
        // Write our fake face files
        std::string output_dir = "HeartGeom3d";
        WriteFakeFaceFile(output_dir, "epi.tri", epi_face, 3u);
        WriteFakeFaceFile(output_dir, "lv.tri", lv_face, 3u);
        WriteFakeFaceFile(output_dir, "rv.tri", rv_face, 3u);

        PetscTools::Barrier(); // Make sure files are written

        // Read in
        OutputFileHandler handler(output_dir, false);
        std::string dir_path = handler.GetOutputDirectoryFullPath();
        HeartGeometryInformation<3> info(mesh, dir_path + "/epi.tri", dir_path + "/lv.tri", dir_path + "/rv.tri", false);

        //covering exceptions
        TS_ASSERT_THROWS_THIS(info.DetermineLayerForEachNode(0.9, 0.9), "The sum of fractions of epicardial and endocardial layers must be lesser than 1");
        TS_ASSERT_THROWS_THIS(info.DetermineLayerForEachNode(0.9, -1.0), "A fraction of a layer must be positive");
        TS_ASSERT_THROWS_THIS(info.DetermineLayerForEachNode(-2.0, 1.0), "A fraction of a layer must be positive");
        
        info.DetermineLayerForEachNode(0.29, 0.51);

        info.WriteLayerForEachNode("TestHeartGeom","layers.het");
        
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            
            if (x < 3 || x > 27)
            {
                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],EPI);
            }
            else if (x < 5 || x > 25)
            {
                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],MID);
            }
            else
            {
                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],ENDO);
            }            
        }

        //now we test the constructor that takes in the node heterogeneity file
        std::string file = OutputFileHandler::GetChasteTestOutputDirectory() + "TestHeartGeom/layers.het";
        HeartGeometryInformation<3> info_2(file);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), info_2.rGetLayerForEachNode().size());
        
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            
            if (x < 3 || x > 27)
            {
                TS_ASSERT_EQUALS(info_2.rGetLayerForEachNode()[index],EPI);
            }
            else if (x < 5 || x > 25)
            {
                TS_ASSERT_EQUALS(info_2.rGetLayerForEachNode()[index],MID);
            }
            else
            {
                TS_ASSERT_EQUALS(info_2.rGetLayerForEachNode()[index],ENDO);
            }   
        }
        
        //covering exceptions
        TS_ASSERT_THROWS_THIS(HeartGeometryInformation<3> info_non_existent_file("rubbish"), "Could not open heterogeneities file (rubbish)");
        TS_ASSERT_THROWS_THIS(HeartGeometryInformation<3> info_bad_file("heart/test/data/ValidPseudoEcg1D.dat"), "A value in the heterogeneities file (heart/test/data/ValidPseudoEcg1D.dat) is out of range (or not an integer). It should be epi = 0, mid = 1, endo = 2");
    }
    
    void TestHeartGeometryTakingFiles() throw (Exception)
    {
        //files containing list of nodes on each surface
        std::string epi_surface = "apps/simulations/propagation3dparallel/heart_chaste2_renum_i_triangles.epi";
        std::string lv_surface = "apps/simulations/propagation3dparallel/heart_chaste2_renum_i_triangles.lv";
        std::string rv_surface = "apps/simulations/propagation3dparallel/heart_chaste2_renum_i_triangles.rv";

               
        //read in the mesh
        TrianglesMeshReader<3,3> mesh_reader("apps/simulations/propagation3dparallel/heart_chaste2_renum_i_triangles");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        //calculate the geometry informations
        HeartGeometryInformation<3> info(mesh, epi_surface, lv_surface, rv_surface, true);
        info.DetermineLayerForEachNode(0.25,0.375);
        //and write them out to file
        OutputFileHandler results_handler("CellularHeterogeneity", false);
        if (PetscTools::AmMaster())
        {        
            out_stream p_file = results_handler.OpenOutputFile("distances.dat");        
            
            for (unsigned index=0; index<mesh.GetNumNodes(); index++)
            {
                (*p_file)<<info.rGetLayerForEachNode()[index]<<std::endl;
            }
            //since we visually checked that the output file is correct,
            //we check that the data in the file match the calculated one.
            EXPECT0(system, "diff " + results_handler.GetOutputDirectoryFullPath() + "/distances.dat " + "heart/test/data/heart_geometry_layers.dat");
        }
    }
};

#endif /*TESTHEARTGEOMETRYINFORMATION_HPP_*/

