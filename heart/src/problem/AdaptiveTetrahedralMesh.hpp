/*

Copyright (C) Fujitsu Laboratories of Europe, 2009-2010

*/

/*

Copyright (C) University of Oxford, 2005-2011

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



#ifndef ADAPTIVETETRAHEDRALMESH_HPP_
#define ADAPTIVETETRAHEDRALMESH_HPP_

#define CXXBLAS_H    // This stops multiple definition of blas headers (one via Chaste, one via libadaptivity/include/cxxblas.h)
                    // Might break things, no idea, will keep it here until problems appear....


#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)

#ifdef CHASTE_ADAPTIVITY
#define HAVE_VTK // Tell libadaptivity that we are using VTK

#include "AbstractTetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"

// Include libadaptivity header files.
// Need to add $LIBADAPTIVITY_DIR/include/ to the include path.
#include "vtk.h"
#include "../metric_field/include/DiscreteGeometryConstraints.h"
#include "../metric_field/include/ErrorMeasure.h"
#include "../adapt3d/include/Adaptivity.h"

/**
 * An adaptive tetrahedral mesh class. Basically just a vtkUnstructuredGrid object with some additional
 * methods and variables to wrap the Imperial College adaptivity library.
 */
class AdaptiveTetrahedralMesh
{
    friend class TestAdaptivityLibrary;
    friend class TestAdaptiveTetrahedralMesh;
    friend class TestAdaptiveTetrahedralMeshLargeMeshes;

private:
    /** vtkUnstructuredGrid: used to interface with adaptivity library */
    vtkUnstructuredGrid *mpVtkUnstructuredGrid;

    /** Number of nodes in the mesh */
    unsigned mNumNodes;

    /** Number of elements in the mesh */
    unsigned mNumElements;

    /** Number of nodes privately owned by this process (i.e. excluding halos) */
    unsigned mNumLocalNodes;

    /** Record whether an adapt has succeeded (i.e. whether or not the adaptivity library has run without returning an error) */
    bool mAdaptSuccess;

    /** DiscreteGeoemtryConstraints object from adaptivity library: identifies co-planar surface elements */
    DiscreteGeometryConstraints* mpDiscreteGeometryConstraints;
    /** ErrorMeasure object from adaptivity library: calculates the metric field */
    ErrorMeasure* mpErrorMeasure;
    /** Adaptivity object from adaptivity library: controls the mesh adaption */
    Adaptivity* mpAdapt;

    /** Adaptivity library requirement */
    std::vector<int> SENList;
    /** Adaptivity library requirement */
    std::vector<int> sids;
    /** Adaptivity library requirement */
    std::vector<double> max_len;

    /**
     * Determine whether or not an edge of the mesh is of sufficient quality. Edge is "good" if the error metric
     * associated  with it is in the range [ 1 - mGoodEdgeRange , 1 + mGoodEdgeRange ]
     */
    double mGoodEdgeRange;

    /** Proportion of edges that must be deemed "bad" (i.e. not good) before an adapt takes place */
    double mBadEdgeCriterion;

    /** Whether or not Adaptivity library should be in verbose mode (default = false) */
    bool mVerbose;

public:

    /**
     * Constructor
     */
    AdaptiveTetrahedralMesh();

    /**
     * Destructor
     */
    ~AdaptiveTetrahedralMesh();

    /**
     * Method to construct an AdaptiveTetrahedralMesh object from a .vtu (vtkUnstructuredGrid format) file
     *
     * @param fileName File name and full path to the file
     */
    void ConstructFromVtuFile(std::string fileName);

    /**
     * Method to construct an AdaptiveTetrahedralMesh object from a Chaste mesh
     *
     * @param rMesh Pointer to the Chaste mesh object
     */
    void ConstructFromMesh(AbstractTetrahedralMesh<3,3>* rMesh);

    /**
     * Method to construct an AdaptiveTetrahedralMesh object from a Chaste distributed mesh
     *
     * @param rMesh Pointer to the Chaste mesh object
     */
    void ConstructFromDistributedMesh(DistributedTetrahedralMesh<3,3>* rMesh);

    /**
     * Add vtkPointData to mpVtkUnstructuredGrid, e.g. to store the values of a variable at each node
     * of the mesh
     *
     * @param dataName Name of the data to be stored
     * @param dataPayload std::vector containing the values at each point
     */
    void AddPointData(std::string dataName, std::vector<double> dataPayload);

    /**
     * Add vtkPointData to mpVtkUnstructuredGrid, e.g. to store the values of a variable at each node
     * of the mesh
     *
     * @param dataName Name of the data to be stored
     * @param dataPayload std::vector containing the values at each point
     */
    void AddPointData(std::string dataName, std::vector<unsigned> dataPayload);

    /**
     * Remove vtkPointData from the mpVtkUnstructuredGrid
     *
     * @param dataName Name of the data to be removed
     */
    void RemoveArray(std::string dataName);

    /**
     * Write out mpVtkUnstructured grid in .vtu format
     *
     * @param directory Directory to write the file in
     * @param fileName File name
     */
    void WriteMeshToFile(std::string directory, std::string fileName);

    /**
     * Write out mpVtkUnstructured grid in .pvtu format
     *
     * @param directory Directory to write the file in
     * @param fileName File name
     */
    void WriteMeshToDistributedFile(std::string directory, std::string fileName);

    /**
     * Return a pointer to mpVtkUnstructuredGrid
     */
    vtkUnstructuredGrid* GetVtkUnstructuredGrid();

    /**
     * Set the values of mGoodEdgeRange and mBadEdgeCriterion to be used in determining whether an adapt is
     * necessary or if the current mesh is of sufficient quality
     *
     * @param range Value of mGoodEdgeRange to be used
     * @param criterion Value of mBadEdgeCriterion to be used
     */
    void SetAdaptCriterion(double range, double criterion);

    /**
     * Return the number of nodes in the mesh
     */
    unsigned GetNumNodes();

    /**
     * Return the number of locally owned nodes, excluding halos
     */
    unsigned GetNumLocalNodes();

    /**
     * Return the number of locally owned nodes, including halos (i.e. mpVtkUnstructuredGrid->GetNumberOfPoints())
     */
    unsigned GetNumLocalAndHaloNodes();

    /**
     * Return the number of elements in the mesh
     */
    unsigned GetNumElements();

    /**
     * Return the number of locally owned elements (i.e. mpVtkUnstructuredGrid->GetNumberOfCells())
     */
    unsigned GetNumLocalElements();

    /**
     * Return the number of surface elements in the mesh. Can only be called after CalculateSENListAndSids(),
     * since vtkUnstructuredGrids do not know about surface elements.
     */
    unsigned GetNumSurfaceElements();

    /**
     * Calculate the surface element-node list and the surface IDs for use in Adaptivity library
     *
     * @param coplanarTolerance Tolerance to be used when determining whether or not two surface elements are coplanar
     */
    void CalculateSENListAndSids(double coplanarTolerance = 0.9999999);

    /**
     * Return the proportion of edges with an error metric value in the range [1 - range, 1 + range].
     *
     * @param range The size of the interval that we are interested in
     */
    double GetEdgeLengthDistribution(double range);

    /**
     * Calculate the discrete geometry constraints for the mesh and an error metric, then adapt the mesh
     * based on these (i.e. without the need to call GetGeometryConstraints() or CalculateErrorMetric())
     *
     * Updates mpVtkUnstructuredGrid to point at a new VTK object that represents the adapted mesh
     *
     */
    void AdaptMesh();

    /**
     * Delete mpDiscreteGeoemtryConstraints, mpErrorMetric and mpAdapt (these need to be created from scratch
     * for each adapt) and clear the entries in SENList, sids and max_len
     */
    void Reset();

    /**
     * Return mAdaptSuccess
     */
    bool GetAdaptSuccess();

    /**
     * Switch to verbose mode
     *
     * @param verbose Bool to specify whether we are switching verbose mode on (true - default) or off (false)
     */
    void MakeVerbose(bool verbose=true);


protected:
    /**
     * Calculate the discrete geometry constraints for the mesh
     *
     * Users should call AdaptMesh() to calculate geometry constraints, error metric and new mesh in one step
     */
    void GetGeometryConstraints();

    /**
     * Calculate an error metric in preparation for adapting (currently uses Vm as the adaptive
     * variable)
     *
     * Users should call AdaptMesh() to calculate geometry constraints, error metric and new mesh in one step
     */
    void CalculateErrorMetric();

    /**
     * Adapt the mesh based on the error metric calculated using CalculateErrorMetric() and subject
     * to the constraints given by GetGeometryConstraints()
     *
     * Updates mpVtkUnstructuredGrid to point at a new VTK object that represents the adapted mesh
     *
     * Users should call AdaptMesh() to calculate geometry constraints, error metric and new mesh in one step
     */
    void Adapt();
};
#endif /*CHASTE_ADAPTIVITY */
#endif /*CHASTE_VTK */

#endif /*ADAPTIVETETRAHEDRALMESH_HPP_*/