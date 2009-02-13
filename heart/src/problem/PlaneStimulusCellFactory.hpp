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


#ifndef PLANESTIMULUSCELLFACTORY_HPP_
#define PLANESTIMULUSCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"
#include "LogFile.hpp"

template<class CELL, unsigned DIM>
class PlaneStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    SimpleStimulus* mpStimulus;

public:
    PlaneStimulusCellFactory(double stimulusMagnitude=-600) : AbstractCardiacCellFactory<DIM>()
    {
        // set the new stimulus
        mpStimulus = new SimpleStimulus(stimulusMagnitude, 0.5);
        LOG(1, "Defined a PlaneStimulusCellFactory<"<<DIM<<"> with SimpleStimulus("<<stimulusMagnitude<<",0.5)\n");
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (this->mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new CELL(this->mpSolver, mpStimulus);
        }
        else
        {
            return new CELL(this->mpSolver, this->mpZeroStimulus);
        }
    }

    ~PlaneStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


#endif /*PLANESTIMULUSCELLFACTORY_HPP_*/
