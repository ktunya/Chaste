#ifndef FACE_HPP_
#define FACE_HPP_

#include "UblasCustomFunctions.hpp"
#include <vector>
#include "Exception.hpp"

template <unsigned DIM>
class Face
{
public:
    /**
     * The vertices of the face, in anticlockwise order. Each vertex must be distinct.
     */
    std::vector< c_vector<double, DIM>* > mVertices;

private:    
    void Increment(typename std::vector< c_vector<double, DIM>* >::iterator& rIterator,
                   Face<DIM>& rFace) const;

public:


    class VertexAndAngle
    {
    public:
        c_vector< double ,DIM >* mpVertex;
        double mAngle;
        bool operator<(const VertexAndAngle& other) const
        {
            return mAngle < other.mAngle;
        }
    };


    /**
     * Compare two faces for equality.
     * Two faces are the same if their vertices differ only by cyclic permutation.
     */
    bool operator==(Face<DIM>& otherFace);
    
    /**
     * Compare two faces for inequality
     */
    bool operator!=(Face<DIM>& otherFace);
    
    /**
     * Return a new face in which the order of the vertices is reversed.
     */
    Face<DIM> operator-();
    
    /**
     * Gets the sum of the length of all edges
     * !!!!! NOTE: Don't use this if you are using a periodic mesh
     * Use GetFacePerimeter(face_index) to takeinto account periodicity.
     */
    double GetPerimeter() const;
    
    /**
     * Gets the area of a voronoi face for 2d space only
     * !!!!! NOTE: Don't use this if you are using a periodic mesh
     * Use GetFaceArea(face_index) to takeinto account periodicity.
     */
    double GetArea() const;
    
    /**
     * Returns number of vertices of a particular face
     */
    unsigned GetNumVertices() const;
    
    /**
     * Returns a vector of vertices
     */
    std::vector< c_vector<double, DIM>* > GetVertices() const;
    
    
    /**
     * Reorders the Vertices of the face anticlockwise
     */
     void OrderVerticesAntiClockwise();
     
     /**
     * @param x x-coordinate
     * @param y y-coordinate
     * @return Polar angle in interval (-PI,PI]
     */
     double ReturnPolarAngle(double x, double y) const;
     
};


#endif /*FACE_HPP_*/
