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


#ifndef FACE_HPP_
#define FACE_HPP_

#include "UblasCustomFunctions.hpp"
#include "Exception.hpp"
#include <vector>


/**
 * A face class for use in the VoronoiTessellation class.
 */
template <unsigned DIM>
class Face
{
private:

    /**
     * Increment the Face vertex iterator.
     * 
     * @param rIterator the Face vertex iterator
     * @param rFace the Face
     */
    void Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                   Face<DIM>& rFace) const;

public:

    /**
     * Helper class containing a pointer to a vertex of the face 
     * and the polar angle from the centre of the face to this 
     * vertex.
     * 
     * \todo This is duplicated in the VoronoiTessellation class; move to a separate file?
     */  
    class VertexAndAngle
    {
    public:

        /** Pointer to a vertex. */
        c_vector<double, DIM>* mpVertex;

        /** Polar angle. */
        double mAngle; 

        /**
         * Less-than angle comparison operator.
         * 
         * @param rOther the VertexAndAngle object to compare to
         */
        bool operator<(const VertexAndAngle& rOther) const
        {
            return mAngle < rOther.mAngle;
        }
    };

    /**
     * The vertices of the face, in anticlockwise order. Each vertex 
     * must be distinct.
     * 
     * This member variable is public as it is accessed directly by 
     * VoronoiTessellation methods.
     */
    std::vector< c_vector<double, DIM>* > mVertices;

    /**
     * Compare two faces for equality. Two faces are the same if their 
     * vertices differ only by cyclic permutation.
     * 
     * @param rOtherFace the Face to compare to
     */
    bool operator==(Face<DIM>& rOtherFace);

    /**
     * Compare two faces for inequality.
     *
     * @param rOtherFace the Face to compare to
     */
    bool operator!=(Face<DIM>& rOtherFace);

    /**
     * Return a new face in which the order of the vertices is reversed.
     */
    Face<DIM> operator-();

    /**
     * Get the sum of the length of all edges of the Face.
     * 
     * NOTE: Don't use this if you are using a periodic mesh. Use 
     * GetFacePerimeter(face_index) to take into account periodicity.
     */
    double GetPerimeter() const;

    /**
     * Get the area of the Face (works in 2D only).
     * 
     * NOTE: Don't use this if you are using a periodic mesh. Use 
     * GetFaceArea(face_index) to takeinto account periodicity.
     */
    double GetArea() const;

    /**
     * Return number of vertices of the Face.
     */
    unsigned GetNumVertices() const;

    /**
     * Return a vector of vertices owned by the Face.
     */
    std::vector< c_vector<double, DIM>* > GetVertices() const;

    /**
     * Reorder the vertices of the Face anticlockwise.
     */
    void OrderVerticesAntiClockwise();

    /**
     * Return the polar angle of the point (x,y).
     * 
     * \todo This is duplicated in the VoronoiTessellation class; move to a separate file?
     * 
     * @param x x-coordinate
     * @param y y-coordinate
     * @return Polar angle in interval (-PI,PI]
     */
    double ReturnPolarAngle(double x, double y) const;

};


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
void Face<DIM>::Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                          Face<DIM>& rFace) const
{
    rIterator++;
    if (rIterator == rFace.mVertices.end())
    {
        rIterator = rFace.mVertices.begin();
    }
}

template<unsigned DIM>
bool Face<DIM>::operator==(Face<DIM>& rOtherFace)
{
    typename std::vector< c_vector<double, DIM>* >::iterator this_iterator = mVertices.begin();
    typename std::vector< c_vector<double, DIM>* >::iterator other_iterator = rOtherFace.mVertices.begin();

    // Find first vertex
    while ( this_iterator!=mVertices.end() &&
            other_iterator!=rOtherFace.mVertices.end() &&
            norm_2(**this_iterator - **other_iterator) >1e-10 )
    {
        this_iterator++;
    }
    if (this_iterator==mVertices.end() || other_iterator==rOtherFace.mVertices.end())
    {
        // Could not find first vertex; faces are distinct unless they are empty
        return ( this_iterator==mVertices.end() && other_iterator==rOtherFace.mVertices.end() );
    }

    typename std::vector< c_vector<double, DIM>* >::iterator this_start=this_iterator;
    Increment(this_iterator, *this);
    Increment(other_iterator, rOtherFace);

    // Check remanining vertices are equal
    while (this_iterator != this_start)
    {
        if (norm_2(**this_iterator - **other_iterator) > 1e-10)
        {
            return false;
        }
        else
        {
            Increment(this_iterator, *this);
            Increment(other_iterator, rOtherFace);
        }
    }
    return (other_iterator == rOtherFace.mVertices.begin());
}

#define COVERAGE_IGNORE // Spuriously not covered
template<unsigned DIM>
bool Face<DIM>::operator!=(Face& rOtherFace)
{
   return !(*this == rOtherFace);
}

#undef COVERAGE_IGNORE

template<unsigned DIM>
Face<DIM> Face<DIM>::operator-()
{
   Face<DIM> reversed_face;
   typename std::vector< c_vector<double, DIM>* >::iterator this_iterator=mVertices.end();
   while (this_iterator != mVertices.begin())
   {
       this_iterator--;
       reversed_face.mVertices.push_back(*this_iterator);
   }
   return reversed_face;
}

template<unsigned DIM>
double Face<DIM>::GetPerimeter() const
{
    double perimeter_return = 0;
    for (unsigned i=0; i<mVertices.size(); i++)
    {
        perimeter_return += norm_2(*mVertices[i]-*mVertices[(i+1)%mVertices.size()]);
    }
    return perimeter_return;
}

template<unsigned DIM>
double Face<DIM>::GetArea() const
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    double area_return = 0;
    for(unsigned i=0; i<mVertices.size(); i++)
    {
        //  Area = sum ( x_i * y_i+1 - y_i * x_i+1 )/2.0 over all vertices,
        //      assuming vertices are ordered anti-clockwise
        area_return += ( (*mVertices[i])(0) * (*mVertices[(i+1)%mVertices.size()])(1)
                        -(*mVertices[i])(1) * (*mVertices[(i+1)%mVertices.size()])(0) ) / 2.0;
    }
    return area_return;
}

template<unsigned DIM>
unsigned Face<DIM>::GetNumVertices() const
{
    return mVertices.size();
}

template<unsigned DIM>
std::vector< c_vector<double, DIM>* > Face<DIM>::GetVertices() const
{
    return mVertices;
}

template<unsigned DIM>
double Face<DIM>::ReturnPolarAngle(double x, double y) const
{
    if (x == 0)
    {
        if (y > 0)
        {
            return M_PI/2.0;
        }
        else if (y < 0)
        {
            return -M_PI/2.0;
        }
        else
        {
            EXCEPTION("Tried to compute polar angle of (0,0)");
        }
    }

    double angle = atan(y/x);

    if (y >= 0 && x < 0 )
    {
        angle += M_PI;
    }
    else if (y < 0 && x < 0 )
    {
        angle -= M_PI;
    }
    return angle;
}

template<unsigned DIM>
void Face<DIM>::OrderVerticesAntiClockwise()
{
    // Reorder mVertices anticlockwise
    std::vector<VertexAndAngle> vertices_and_angles;

    c_vector<double,DIM> centre = zero_vector<double>(DIM);

    for (unsigned j=0; j<mVertices.size(); j++)
    {
        centre += *(mVertices[j]);
    }

    centre /= mVertices.size();
    for (unsigned j=0; j<mVertices.size(); j++)
    {
        VertexAndAngle va;
        c_vector<double, DIM> centre_to_vertex = *(mVertices[j]) - centre;

        va.mAngle = ReturnPolarAngle(centre_to_vertex(0), centre_to_vertex(1));
        va.mpVertex = mVertices[j];
        vertices_and_angles.push_back(va);
    }

    std::sort(vertices_and_angles.begin(), vertices_and_angles.end());

    // Create face
    mVertices.clear();
    for (typename std::vector<VertexAndAngle>::iterator vertex_iterator = vertices_and_angles.begin();
         vertex_iterator !=vertices_and_angles.end();
         vertex_iterator++)
    {
        mVertices.push_back(vertex_iterator->mpVertex);
    }
}

#endif /*FACE_HPP_*/
