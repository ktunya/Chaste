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


#include "SolidMechanicsProblemDefinition.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"


template<unsigned DIM>
const double SolidMechanicsProblemDefinition<DIM>::FREE = DBL_MAX;

template<unsigned DIM>
SolidMechanicsProblemDefinition<DIM>::SolidMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh)
    : mrMesh(rMesh),
      mDensity(1.0),
      mBodyForceType(CONSTANT_BODY_FORCE),
      mConstantBodyForce(zero_vector<double>(DIM)),
      mTractionBoundaryConditionType(NO_TRACTIONS)
{
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetDensity(double density)
{
    assert(density>0.0);
    mDensity = density;
}

template<unsigned DIM>
double SolidMechanicsProblemDefinition<DIM>::GetDensity()
{
    return mDensity;
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetBodyForce(c_vector<double,DIM> bodyForce)
{
    mBodyForceType = CONSTANT_BODY_FORCE;
    mConstantBodyForce = bodyForce;
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetBodyForce(c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t))
{
    mBodyForceType = FUNCTIONAL_BODY_FORCE;
    mpBodyForceFunction = pFunction;
}

template<unsigned DIM>
BodyForceType SolidMechanicsProblemDefinition<DIM>::GetBodyForceType()
{
    return mBodyForceType;
}

template<unsigned DIM>
c_vector<double,DIM> SolidMechanicsProblemDefinition<DIM>::GetConstantBodyForce()
{
    assert(mBodyForceType==CONSTANT_BODY_FORCE);
    return mConstantBodyForce;
}

template<unsigned DIM>
c_vector<double,DIM> SolidMechanicsProblemDefinition<DIM>::EvaluateBodyForceFunction(c_vector<double,DIM>& rX, double t)
{
    assert(mBodyForceType==FUNCTIONAL_BODY_FORCE);
    return (*mpBodyForceFunction)(rX,t);
}


template<unsigned DIM>
TractionBoundaryConditionType SolidMechanicsProblemDefinition<DIM>::GetTractionBoundaryConditionType()
{
    return mTractionBoundaryConditionType;
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                                                         std::vector<c_vector<double,DIM> >& rElementwiseTractions)
{

    assert(rTractionBoundaryElements.size()==rElementwiseTractions.size());
    mTractionBoundaryConditionType = ELEMENTWISE_TRACTION;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mElementwiseTractions = rElementwiseTractions;
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*> rTractionBoundaryElements,
                                                                         c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t))
{
    mTractionBoundaryConditionType=FUNCTIONAL_TRACTION;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mpTractionBoundaryConditionFunction = pFunction;
}


template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetApplyNormalPressureOnDeformedSurface(std::vector<BoundaryElement<DIM-1,DIM>*> rTractionBoundaryElements,
                                                                                   double normalPressure)
{
    mTractionBoundaryConditionType = PRESSURE_ON_DEFORMED;
    mTractionBoundaryElements = rTractionBoundaryElements;
    mNormalPressure = normalPressure;

}



template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetZeroDisplacementNodes(std::vector<unsigned>& rFixedNodes)
{
    mFixedNodes = rFixedNodes;

    for (unsigned i=0; i<mFixedNodes.size(); i++)
    {
        assert(mFixedNodes[i] < mrMesh.GetNumNodes());
    }

    mFixedNodeDisplacements.clear();
    for (unsigned i=0; i<mFixedNodes.size(); i++)
    {
        mFixedNodeDisplacements.push_back(zero_vector<double>(DIM));
    }
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetFixedNodes(std::vector<unsigned>& rFixedNodes, std::vector<c_vector<double,DIM> >& rFixedNodeLocations)
{
    assert(rFixedNodes.size()==rFixedNodeLocations.size());
    mFixedNodes = rFixedNodes;

    mFixedNodeDisplacements.clear();
    for (unsigned i=0; i<mFixedNodes.size(); i++)
    {
        unsigned index = mFixedNodes[i];
        c_vector<double,DIM> displacement;
        for(unsigned j=0; j<DIM; j++)
        {
        	double location = rFixedNodeLocations[i](j);

			// compute the displacement, assuming the node
			// is not free in this direction
        	if(location != FREE)
        	{
        		displacement(j) = location - mrMesh.GetNode(index)->rGetLocation()[j];
            }
            else
            {
            	displacement(j) = FREE;
            }
        }
    	mFixedNodeDisplacements.push_back(displacement);
    }
}


template<unsigned DIM>
std::vector<unsigned>& SolidMechanicsProblemDefinition<DIM>::rGetFixedNodes()
{
    return mFixedNodes;
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& SolidMechanicsProblemDefinition<DIM>::rGetFixedNodeDisplacements()
{
    return mFixedNodeDisplacements;
}

template<unsigned DIM>
std::vector<BoundaryElement<DIM-1,DIM>*>& SolidMechanicsProblemDefinition<DIM>::rGetTractionBoundaryElements()
{
    return mTractionBoundaryElements;
}


template<unsigned DIM>
std::vector<c_vector<double,DIM> >& SolidMechanicsProblemDefinition<DIM>::rGetElementwiseTractions()
{
    assert(mTractionBoundaryConditionType==ELEMENTWISE_TRACTION);
    return mElementwiseTractions;
}


template<unsigned DIM>
double SolidMechanicsProblemDefinition<DIM>::GetNormalPressure()
{
    assert(mTractionBoundaryConditionType==PRESSURE_ON_DEFORMED);
    return mNormalPressure;
}

template<unsigned DIM>
c_vector<double,DIM> SolidMechanicsProblemDefinition<DIM>::EvaluateTractionFunction(c_vector<double,DIM>& rX, double t)
{
    assert(mTractionBoundaryConditionType==FUNCTIONAL_TRACTION);
    return (*mpTractionBoundaryConditionFunction)(rX,t);
}



template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetHomogenenousMaterial(AbstractMaterialLaw<DIM>* pMaterialLaw)
{
    mIsHeterogeneousMaterial = false;
    mMaterialLaws.clear();
    mMaterialLaws.push_back(pMaterialLaw);
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::SetHeterogeneousMaterial(std::vector<AbstractMaterialLaw<DIM>*>& rMaterialLaws)
{
    mIsHeterogeneousMaterial = true;
    mMaterialLaws.clear();
    assert(mrMesh.GetNumElements()==rMaterialLaws.size());
    mMaterialLaws = rMaterialLaws;
}

template<unsigned DIM>
bool SolidMechanicsProblemDefinition<DIM>::IsHeterogeneousMaterial()
{
    assert(mMaterialLaws.size()!=0);
    return mIsHeterogeneousMaterial;
}

template<unsigned DIM>
AbstractMaterialLaw<DIM>* SolidMechanicsProblemDefinition<DIM>::GetLawForHomogeneousMaterial()
{
    assert(!mIsHeterogeneousMaterial);
    assert(mMaterialLaws.size()==1);
    return mMaterialLaws[0];
}

template<unsigned DIM>
AbstractMaterialLaw<DIM>* SolidMechanicsProblemDefinition<DIM>::GetLawForHeterogeneousMaterial(unsigned elementIndex)
{
    assert(mIsHeterogeneousMaterial);
    assert(mMaterialLaws.size()==mrMesh.GetNumElements());
    assert(elementIndex < mrMesh.GetNumElements());
    return mMaterialLaws[elementIndex];
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::VerifyIncompressibleMaterialLaws()
{
    if(mMaterialLaws.size()==0)
    {
        EXCEPTION("No material laws have been set on SolidMechanicsProblemDefinition");
    }

    for (unsigned i=0; i<mMaterialLaws.size(); i++)
    {
        AbstractIncompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractIncompressibleMaterialLaw<DIM>*>(mMaterialLaws[i]);
        if (!p_law)
        {
            EXCEPTION("SolidMechanicsProblemDefinition used in incompressible problem but has been given compressible material law(s)");
        }
    }
}

template<unsigned DIM>
void SolidMechanicsProblemDefinition<DIM>::VerifyCompressibleMaterialLaws()
{
    if(mMaterialLaws.size()==0)
    {
        EXCEPTION("No material laws have been set on SolidMechanicsProblemDefinition");
    }

    for (unsigned i=0; i<mMaterialLaws.size(); i++)
    {
        AbstractCompressibleMaterialLaw<DIM>* p_law = dynamic_cast<AbstractCompressibleMaterialLaw<DIM>*>(mMaterialLaws[i]);
        if (!p_law)
        {
            EXCEPTION("SolidMechanicsProblemDefinition used in compressible problem but has been given incompressible material law(s)");
        }
    }
}


//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class SolidMechanicsProblemDefinition<2>;
template class SolidMechanicsProblemDefinition<3>;
