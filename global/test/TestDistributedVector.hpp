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


#ifndef TESTDISTRIBUTEDVECTOR_HPP_
#define TESTDISTRIBUTEDVECTOR_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

#include "DistributedVector.hpp"
#include "PetscTools.hpp"

#include "DistributedVectorFactory.hpp"

#include "OutputFileHandler.hpp"

class TestDistributedVector : public CxxTest::TestSuite
{
public:
    void TestDistributedVectorFactory()
    {
        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned total=100;
        DistributedVectorFactory factory(total);
        PetscInt petsc_lo, petsc_hi;
        unsigned expected_local_size = (total/num_procs);
        unsigned max_expected_local_size = (total/num_procs) + total%num_procs;
         
        Vec petsc_vec;
        petsc_vec = factory.CreateVec();
        VecGetOwnershipRange(petsc_vec, &petsc_lo, &petsc_hi);
        unsigned local_size = petsc_hi - petsc_lo;
        if (expected_local_size != local_size)
        {
            //This test is not robust.  PETSc has been known to split a 100 elements among 7 processors like this
            // 0  1  2  3  4  5  6  7
            //15 15 14 14 14 14 14 14 
            //TS_ASSERT_EQUALS(local_size, max_expected_local_size);
            //This test will fail if there are fewer elements in the vector than the number of processors
            TS_ASSERT_LESS_THAN_EQUALS(local_size, max_expected_local_size);
        }    
        //Test that we can construct a factory from a given vector
        DistributedVectorFactory factory2(petsc_vec);
        Vec petsc_vec2;
        petsc_vec2 = factory2.CreateVec();
        
        VecGetOwnershipRange(petsc_vec2, &petsc_lo, &petsc_hi);
        unsigned local_size2 = petsc_hi - petsc_lo;
        TS_ASSERT_EQUALS(local_size, local_size2);
        TS_ASSERT_EQUALS(local_size, factory2.GetLocalOwnership());
        TS_ASSERT_EQUALS(total, factory2.GetProblemSize());
        
        TS_ASSERT_EQUALS((unsigned)(petsc_hi), factory2.GetHigh());
        TS_ASSERT_EQUALS((unsigned)(petsc_lo), factory2.GetLow());
        

        //Uneven test (as above).
        //Calculate total number of elements in the vector
        unsigned total_elements = (num_procs+1)*num_procs/2;
        unsigned my_rank=PetscTools::GetMyRank();

        DistributedVectorFactory uneven_factory(total_elements, my_rank+1);

        Vec petsc_vec_uneven = uneven_factory.CreateVec();

        VecGetOwnershipRange(petsc_vec_uneven, &petsc_lo, &petsc_hi);

        unsigned expected_lo = (my_rank+1)*my_rank/2;
        unsigned expected_hi = (my_rank+2)*(my_rank+1)/2;

        TS_ASSERT_EQUALS((unsigned)petsc_lo, expected_lo);
        TS_ASSERT_EQUALS((unsigned)petsc_hi, expected_hi);
        
        VecDestroy(petsc_vec);
        VecDestroy(petsc_vec2);
        VecDestroy(petsc_vec_uneven);  
    }

    void TestRead()
    {
        // WRITE VECTOR
        // create a 10 element petsc vector
        unsigned vec_size = 10u;
        Vec vec;
        VecCreate(PETSC_COMM_WORLD, &vec);
        VecSetSizes(vec, PETSC_DECIDE, vec_size);
        VecSetFromOptions(vec);
        // calculate the range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(vec, &petsc_lo, &petsc_hi);
        unsigned lo=(unsigned)petsc_lo;
        unsigned hi=(unsigned)petsc_hi;
        // create 20 element petsc vector
        Vec striped;
        VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*vec_size, &striped);
        // write some values
        double* p_vec;
        VecGetArray(vec, &p_vec);
        double* p_striped;
        VecGetArray(striped, &p_striped);
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            p_vec[local_index] = local_index*global_index;
            p_striped[2*local_index  ] = local_index;
            p_striped[2*local_index+1] = global_index*global_index;
        }
        VecRestoreArray(vec, &p_vec);
        VecAssemblyBegin(vec);
        VecAssemblyEnd(vec);
        VecRestoreArray(striped, &p_striped);
        VecAssemblyBegin(striped);
        VecAssemblyEnd(striped);

        // READ VECTOR
        DistributedVectorFactory factory(vec);
        DistributedVector distributed_vector = factory.CreateDistributedVector(vec);
        DistributedVector distributed_vector2 = factory.CreateDistributedVector(striped);
        DistributedVector::Stripe linear(distributed_vector2,0);
        DistributedVector::Stripe quadratic(distributed_vector2,1);
        // check the range
        TS_ASSERT_EQUALS(factory.GetProblemSize(), vec_size);
        TS_ASSERT_EQUALS(distributed_vector.Begin().Global,lo);
        TS_ASSERT_EQUALS(distributed_vector.End().Global,hi);
        // read some values
        for (DistributedVector::Iterator index = distributed_vector.Begin();
             index!= distributed_vector.End();
             ++index)
        {
            TS_ASSERT_EQUALS(distributed_vector[index], index.Local*index.Global);
            TS_ASSERT_EQUALS(linear[index], index.Local);
            TS_ASSERT_EQUALS(quadratic[index], index.Global * index.Global);
        }

        // read the 2nd element of the first vector
        if (lo<=2 && 2<hi)
        {
            TS_ASSERT(distributed_vector.IsGlobalIndexLocal(2));
            TS_ASSERT_EQUALS(distributed_vector[2],2*(2-lo));
        }
        else
        {
            TS_ASSERT(!distributed_vector.IsGlobalIndexLocal(2));
            TS_ASSERT_THROWS_ANYTHING(distributed_vector[2]);
        }

        //read the 3rd element of the other vectors
        if (lo<=3 && 3<hi)
        {
            TS_ASSERT(distributed_vector.IsGlobalIndexLocal(3));
            TS_ASSERT_EQUALS(linear[3],(3-lo));
            TS_ASSERT_EQUALS(quadratic[3],3*3);
        }
        else
        {
            TS_ASSERT(!distributed_vector.IsGlobalIndexLocal(3));
            TS_ASSERT_THROWS_ANYTHING(linear[3]);
            TS_ASSERT_THROWS_ANYTHING(quadratic[3]);
        }

        VecDestroy(vec);
        VecDestroy(striped);
    }

    void TestWrite()
    {
        //WRITE VECTOR
        // create a 10 element petsc vector
        DistributedVectorFactory factory(10);
        Vec striped=factory.CreateVec(2);
        Vec chunked=factory.CreateVec(2);
        Vec petsc_vec=factory.CreateVec();

        DistributedVector distributed_vector = factory.CreateDistributedVector(petsc_vec);
        DistributedVector distributed_vector_striped = factory.CreateDistributedVector(striped);
        DistributedVector distributed_vector_chunked = factory.CreateDistributedVector(chunked);
        DistributedVector::Stripe linear(distributed_vector_striped, 0);
        DistributedVector::Stripe quadratic(distributed_vector_striped, 1);
        DistributedVector::Chunk linear_chunk(distributed_vector_chunked, 0);
        DistributedVector::Chunk quadratic_chunk(distributed_vector_chunked, 1);
        // write some values
        for (DistributedVector::Iterator index = distributed_vector.Begin();
             index!= distributed_vector.End();
             ++index)
        {
            distributed_vector[index] =  -(double)(index.Local*index.Global);
            linear[index] =  -1;
            quadratic[index] =  index.Local+1;
            linear_chunk[index] = -1;
            quadratic_chunk[index] =  index.Global+1;
        }

        distributed_vector.Restore();
        distributed_vector_striped.Restore();
        distributed_vector_chunked.Restore();

        //READ VECTOR
        // calculate my range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(petsc_vec,&petsc_lo,&petsc_hi);
        unsigned lo=(unsigned)petsc_lo;
        unsigned hi=(unsigned)petsc_hi;
        // read some values
        double* p_striped;
        VecGetArray(striped, &p_striped);
        double* p_chunked;
        VecGetArray(chunked, &p_chunked);
        double* p_vec;
        VecGetArray(petsc_vec, &p_vec);
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            TS_ASSERT_EQUALS(p_vec[local_index], -(double)local_index*global_index);
            TS_ASSERT_EQUALS(p_striped[2*local_index], -1.0);
            TS_ASSERT_EQUALS(p_striped[2*local_index+1], local_index+1);

            TS_ASSERT_EQUALS(linear[global_index], -1.0);
            TS_ASSERT_EQUALS(quadratic[global_index], local_index+1);

            TS_ASSERT_EQUALS(p_chunked[local_index], -1.0);
            TS_ASSERT_EQUALS(p_chunked[ (hi - lo) + local_index], global_index+1);

            TS_ASSERT_EQUALS(linear_chunk[global_index], -1.0);
            TS_ASSERT_EQUALS(quadratic_chunk[global_index], global_index+1);
        }

        //Read item 2 from the distributed vectors (for coverage)
        if (lo<=2 && 2<hi)
        {
            TS_ASSERT(distributed_vector.IsGlobalIndexLocal(2));
            TS_ASSERT_EQUALS(linear[2], -1.0);
            TS_ASSERT_EQUALS(quadratic[2], 3.0 - lo);
            TS_ASSERT_EQUALS(linear_chunk[2], -1.0);
            TS_ASSERT_EQUALS(quadratic_chunk[2], 3.0);
        }
        else
        {
            TS_ASSERT(!distributed_vector.IsGlobalIndexLocal(2));
            TS_ASSERT_THROWS_ANYTHING(linear[2]);
            TS_ASSERT_THROWS_ANYTHING(linear_chunk[2]);
        }

        VecDestroy(petsc_vec);
        VecDestroy(striped);
        VecDestroy(chunked);
    }

    void TestException()
    {
        TS_ASSERT_THROWS_ANYTHING(throw DistributedVectorException() );
    }

    void TestUnevenDistribution()
    {
        unsigned my_rank = PetscTools::GetMyRank();
        unsigned num_procs = PetscTools::GetNumProcs();

        //Calculate total number of elements in the vector
        unsigned total_elements = (num_procs+1)*num_procs/2;

        DistributedVectorFactory factory(total_elements, my_rank+1);
        Vec petsc_vec = factory.CreateVec(1);

        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(petsc_vec,&petsc_lo,&petsc_hi);

        unsigned expected_lo = (my_rank+1)*my_rank/2;
        unsigned expected_hi = (my_rank+2)*(my_rank+1)/2;

        TS_ASSERT_EQUALS((unsigned)petsc_lo, expected_lo);
        TS_ASSERT_EQUALS((unsigned)petsc_hi, expected_hi);

        VecDestroy(petsc_vec);
    }
    
    void TestArchiving()
    {
        const unsigned TOTAL = 100;
        DistributedVectorFactory factory(TOTAL);
        unsigned num_local_items = factory.GetLocalOwnership();
        unsigned hi = factory.GetHigh();
        unsigned lo = factory.GetLow();
        TS_ASSERT_EQUALS(factory.GetProblemSize(), TOTAL);
        
        // Where to archive
        OutputFileHandler handler("archive");
        std::string archive_filename = handler.GetProcessUniqueFilePath("factory.arch");
        
        // Archive
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            const DistributedVectorFactory* const p_factory = &factory;
            output_arch << p_factory;
        }
        
        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            DistributedVectorFactory* p_new_factory;
            input_arch >> p_new_factory;
            
            TS_ASSERT_EQUALS(p_new_factory->GetProblemSize(), TOTAL);
            TS_ASSERT_EQUALS(p_new_factory->GetHigh(), hi);
            TS_ASSERT_EQUALS(p_new_factory->GetLow(), lo);
            TS_ASSERT_EQUALS(p_new_factory->GetLocalOwnership(), num_local_items);
            
            delete p_new_factory;
        }
        
        //Restore from a single process archive
        {
            std::ifstream ifs("global/test/data/distributed_vector_factory.arch", std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            DistributedVectorFactory* p_new_factory = NULL;
            
            if (PetscTools::IsSequential())
            {
                input_arch >> p_new_factory;
            
                TS_ASSERT_EQUALS(p_new_factory->GetProblemSize(), TOTAL);
                TS_ASSERT_EQUALS(p_new_factory->GetHigh(), TOTAL);
                TS_ASSERT_EQUALS(p_new_factory->GetLow(), 0U);
                TS_ASSERT_EQUALS(p_new_factory->GetLocalOwnership(), TOTAL);
            }
            else
            {
                //Should not read this archive
                TS_ASSERT_THROWS_ANYTHING(input_arch >> p_new_factory);
            }
            
            delete p_new_factory;
        }
    }
};

#endif /*TESTDISTRIBUTEDVECTOR_HPP_*/
