#!/usr/bin/env python

"""Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
 
"""
Convert a .ele .cable .edge or .face file which is in binary format and has (integer) attributes

This is needed because 
 * Chaste 3.2 and later use 8-byte double precision floats for element attributes
 * Chaste 3.1 and previous used 4-byte unsigned integers for element attributes

Note that Ascii files are unaffected by the change because the integer is automatically converted on reading

"""
import fileinput
import sys
import struct
import os

if __name__ == "__main__":
    # Checking command line arguments
    if len(sys.argv) != 3:
        print >> sys.stderr, "Usage:", sys.argv[0], "<input_mesh_file> <output_mesh_file>"
        sys.exit(1)
    # Reading in command line arguments
    input_name = sys.argv[1]
    output_name = sys.argv[2]
    ext = os.path.splitext(input_name)[1]

    if not (ext == '.ele' or ext == '.cable' or ext == '.edge' or ext == '.face' ):
        print >> sys.stderr, "This utility is only needed on .ele .cable .edge and .face files which are in binary format"
        sys.exit(2)
    in_file = open(input_name, 'rb')
    head_line =  in_file.readline()
    headers = head_line.split()
    if headers[-1] != 'BIN':
        print >> sys.stderr, "This utility is only needed on .ele .cable .edge and .face files which are in binary format"
        sys.exit(3)
    assert( len(headers) == 4)
    num_eles = int(headers[0])
    nodes_per_element = int(headers[1])
    num_attributes = int(headers[2])
    assert( num_attributes == 1) # Todo - exit gracefully if there's nothing to do?
    old_data_width = nodes_per_element*4 + num_attributes*4
    old_data_size = old_data_width * num_eles + 2
    new_data_width = nodes_per_element*4 + num_attributes*8
    new_data_size = new_data_width * num_eles + 2
    
    # Remember the place on disk
    data_start = in_file.tell()
    # TODO: This doesn't work because \n == 0a == 10 is almost certain to appear in the data
    size_on_disk = len(in_file.readline())  # Todo Does this use memory?  If so, we may as well read once and keep it in memory
    #print size_on_disk, old_data_size, new_data_size
    if size_on_disk == new_data_size:
        print >> sys.stderr, "This file looks like it has already been converted"
        sys.exit(4)
    assert( size_on_disk == old_data_size )
    
    # Write out the header line (which includes the \n)
    out_file = open(output_name, 'w')
    print >>out_file, head_line,

    # Read and convert an element at a time
    in_file.seek(data_start)
    for i in range(0, num_eles):
        # Copy across the element's node indices
        out_file.write( in_file.read(nodes_per_element*4) )
        attribute = struct.unpack('I', in_file.read(4))[0] # Read as a uint 'I'
        assert(attribute == 0) # I'd like to monitor this!
        # Tweak it
        out_file.write(struct.pack('d', attribute)) # Written as double 'd'
    # Repeat closing #\n
    out_file.write(in_file.readline())
    # Repeat provenance
    out_file.write(in_file.readline().strip('\n'))
    # Give a message
    print >>out_file, ' **Attributes converted to Chaste 3.2 form**\n'
    
        

    






