"""Host-specific configuration.

This module contains the logic to pick up the right configuration file
for the machine it is run on, based on the machine's hostname.

Each configuration file should be a Python program providing certain variables
specifying where to find libraries and tools.  These variables are:
 * petsc_2_2_path - path to PETSc 2.2 install (None if not present)
 * petsc_2_3_path - path to PETSc 2.3 install (None if not present)
 * dealii_path    - path to Deal.II install (None if not present)
 * metis_path     - path to METIS install (None if not present)
 * intel_path     - path to Intel compiler installation

 * other_includepaths - list of paths containing other header files
 * other_libpaths     - list of paths containing other libraries, including metis, xsd, and boost
 * other_libraries    - list of other libraries to link against.  This *must* include
                        lapack, blas, and boost_serialization, as their names vary across systems.

 * ccflags - any extra compiler flags needed, as a string

 * tools - a dictionary mapping executable names to their absolute paths, for tools
           not found on $PATH.

Any non-absolute paths will be considered relative to the root of the Chaste install.
"""

import os
import socket
import sys

machine_fqdn = socket.getfqdn()


if machine_fqdn in ["userpc30.comlab.ox.ac.uk", "userpc33.comlab.ox.ac.uk"]:
    import joe as conf
elif machine_fqdn in ["userpc44.comlab.ox.ac.uk", "userpc60.comlab.ox.ac.uk",
                      "userpc58.comlab.ox.ac.uk", "userpc59.comlab.ox.ac.uk"]:
    import chaste as conf
elif machine_fqdn == "zuse.osc.ox.ac.uk":
    import zuse as conf
elif machine_fqdn.endswith(".comlab.ox.ac.uk"):
    import comlab as conf
elif machine_fqdn.startswith('finarfin'):
    import finarfin as conf
elif machine_fqdn.endswith(".maths.nottingham.ac.uk"):
    import nottingham as conf
elif machine_fqdn.endswith(".maths.ox.ac.uk"):
    import maths as conf
else:
    print >>sys.stderr, "Unrecognised machine %s; please add a stanza for it to hostconfig.py" % machine_fqdn
    sys.exit(1)

# For debugging
#for name in dir(conf):
#    if name[0] != '_':
#        print name, '=', getattr(conf, name)

# This is a bit ugly at present: SConstruct calls configure() to fill
# these global variables in, then reads them directly.
# Note that the order in which things are added to these lists often matters!
libpaths = []
incpaths = []
libraries = []

def do_petsc(version, optimised, profile=False, production=False):
    """Determine PETSc include and library paths.

    The locations vary depending on the version of PETSc, and possibly
    whether optimised libraries are to be used.

    The version can be given as 2_2 or 2_3 to choose PETSc minor version.
    If a host doesn't support 2.3, we attempt to use 2.2 instead.  A
    ValueError is raised if 2.2 isn't present when asked for.

    Set optimised to True to use optimised builds of the libraries rather
    than debug builds.
    Set profile to True to use profile builds of PETSc.
    """
    if version == '2_3' and conf.petsc_2_3_path is None:
        # Use 2.2 instead
        version = '2_2'
    if version == '2_2' and conf.petsc_2_2_path is None:
        # Raise a friendly error
        raise ValueError('PETSc 2.2 required, but no path given in the host config.')
    if version == '2_2':
        petsc_base = os.path.abspath(conf.petsc_2_2_path)
        # Gracefully fall back to optimised/non-opt if the requested one isn't there
        if optimised:
            dirs = ['libO_c++', 'libg_c++']
        else:
            dirs = ['libg_c++', 'libO_c++']
        for d in dirs:
            libpath = os.path.join(petsc_base, 'lib', d, conf.petsc_build_name)
            if os.path.exists(libpath): break
        else:
            raise ValueError('No PETSc 2.2 libraries found.')
        incpaths.append(os.path.join(petsc_base, 'bmake', conf.petsc_build_name))
    else:
        petsc_base = os.path.abspath(conf.petsc_2_3_path)
        if production:
            build_name = conf.petsc_build_name_production
        elif profile:
            optimised = False
            build_name = conf.petsc_build_name_profile
        elif optimised:
            build_name = conf.petsc_build_name_optimized
        else:
            build_name = conf.petsc_build_name
        libpath = os.path.join(petsc_base, 'lib', build_name)
        incpaths.append(os.path.join(petsc_base, 'bmake', build_name))
    incpaths.append(os.path.join(petsc_base, 'include'))
    libpaths.append(libpath)
    libraries.extend(['petscts', 'petscsnes', 'petscksp', 'petscdm', 
                      'petscmat', 'petscvec', 'petsc'])

def do_metis():
    """Add METIS include and library paths."""
    if conf.metis_path is None:
        raise ValueError('METIS required, but no path given in the host config.')
    libpath = os.path.abspath(conf.metis_path)
    incpath = os.path.join(libpath, 'Lib') # Yes, METIS is odd!
    libpaths.append(libpath)
    incpaths.append(incpath)
    libraries.append('metis')

def do_dealii(build):
    """Add Deal.II include & library paths, and libraries.

    Deal.II uses different library *names* to distinguish optimised versions.
    """
    if conf.dealii_path is None:
        raise ValueError('Deal.II required, but no path given in the host config.')
    base = os.path.abspath(conf.dealii_path)
    libpaths.append(os.path.join(base, 'lib'))
    relative_incpaths = ['base/include', 'lac/include', 'deal.II/include']
    incpaths.extend(map(lambda relpath: os.path.join(base, relpath),
                        relative_incpaths))
    libs = ['deal_II_1d', 'deal_II_2d', 'deal_II_3d', 'lac', 'base']
    if build.dealii_debugging:
        libs = map(lambda s: s + '.g', libs)
    libraries.extend(libs)

def configure(build):
    """Given a build object (BuildTypes.BuildType instance), configure the build."""
    if build.using_dealii:
        do_petsc('2_2', build.is_optimised) # Deal.II only supports PETSc 2.2
        do_dealii(build)
        do_metis()
        libraries.extend(['blas', 'lapack']) # Use versions provided with Deal.II
    else:
        do_petsc('2_3', build.is_optimised, build.is_profile, build.is_production)
        if build.is_production:
            libraries.extend(conf.blas_lapack_production)
        else:
            libraries.extend(conf.blas_lapack)
    if build.CompilerType() == 'intel':
        intel_path = os.path.abspath(conf.intel_path)
        libpaths.append(os.path.join(intel_path, 'lib'))
    incpaths.extend(conf.other_includepaths)
    libpaths.extend(map(os.path.abspath, conf.other_libpaths))
    libraries.extend(conf.other_libraries)

    build.tools.update(conf.tools)

    if build.CompilerType() == 'intel':
        # Switch to use Intel toolchain
        build.tools['mpicxx'] += ' -CC=icpc'
        build.tools['cxx'] = os.path.join(intel_path, 'bin', 'icpc')
        build.tools['ar'] = os.path.join(intel_path, 'bin', 'xiar')

def ccflags():
    try:
        return conf.ccflags
    except:
        return ''
