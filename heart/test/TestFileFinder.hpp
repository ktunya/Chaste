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

#ifndef TESTFILEFINDER_HPP_
#define TESTFILEFINDER_HPP_

#include <cxxtest/TestSuite.h>
#include "FileFinder.hpp"
#include "ChasteBuildRoot.hpp"
#include "OutputFileHandler.hpp"

class TestFileFinder : public CxxTest::TestSuite
{
public:
    void TestFileFinderOpening() throw(Exception)
    {
        {
            // Can we find our own source file?
            std::string file_name = "heart/src/io/FileFinder.hpp";
            cp::path_type path(file_name);
            path.relative_to(cp::relative_to_type::chaste_source_root);
            FileFinder file_finder(path);
            TS_ASSERT(file_finder.Exists());
            // Check the path is as expected
            std::string abs_path = ChasteBuildRootDir() + file_name;
            TS_ASSERT_EQUALS(file_finder.GetAbsolutePath(), abs_path);
            
            // CWD should be the Chaste source root
            path.relative_to(cp::relative_to_type::cwd);
            FileFinder file_finder2(path);
            TS_ASSERT(file_finder2.Exists());
            // Check the path is as expected
            TS_ASSERT_EQUALS(file_finder2.GetAbsolutePath(), abs_path);
            
            // Alternative constructor
            FileFinder file_finder3(file_name, cp::relative_to_type::chaste_source_root);
            TS_ASSERT(file_finder3.Exists());
            TS_ASSERT_EQUALS(file_finder3.GetAbsolutePath(), abs_path);
        }
        
        {
            // Now check a file in the output directory
            std::string dir_name = "TestFileFinder";
            OutputFileHandler handler(dir_name);
            std::string file_name = "TestFile";
            cp::path_type path(dir_name + "/" + file_name);
            path.relative_to(cp::relative_to_type::chaste_test_output);
            FileFinder file_finder(path);
            TS_ASSERT(! file_finder.Exists());
            // Check the path is as expected
            std::string abs_path = handler.GetOutputDirectoryFullPath() + file_name;
            TS_ASSERT_EQUALS(file_finder.GetAbsolutePath(), abs_path);
            // Create the file
            out_stream fp = handler.OpenOutputFile(file_name);
            fp->close();
            TS_ASSERT(file_finder.Exists());
            
            // Check when providing an absolute path
            cp::path_type path_abs(abs_path);
            path_abs.relative_to(cp::relative_to_type::absolute);
            FileFinder file_finder2(path_abs);
            TS_ASSERT(file_finder2.Exists());
            TS_ASSERT_EQUALS(file_finder2.GetAbsolutePath(), abs_path);
            
            // Alternative constructor
            FileFinder file_finder3(abs_path, cp::relative_to_type::absolute);
            TS_ASSERT(file_finder3.Exists());
            TS_ASSERT_EQUALS(file_finder3.GetAbsolutePath(), abs_path);
        }
    }
    
    void TestNewer()
    {
        FileFinder file("heart/src/io/FileFinder.hpp", cp::relative_to_type::chaste_source_root);
        // A file can't be newer than itself
        TS_ASSERT(!file.IsNewerThan(file));
        // A newly created file better be newer than ourself!
        OutputFileHandler handler("TestFileFinder");
        out_stream fp = handler.OpenOutputFile("new_file");
        fp->close();
        FileFinder new_file("TestFileFinder/new_file", cp::relative_to_type::chaste_test_output);
        TS_ASSERT(new_file.IsNewerThan(file));
        TS_ASSERT(!file.IsNewerThan(new_file));
    }
};

#endif /*TESTFILEFINDER_HPP_*/


