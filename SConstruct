#
# SConstruct
# 
import os

home_dir = os.getcwd()


# CPPPATH: sets the include path for C++
env = Environment()
env.Append(CPPPATH=['/vol/software/software/tools/tmv/tmv0.72/x86_64/include'])
env.Append(CPPPATH=['/opt/local/include'])
env.Append(CPPPATH=[os.path.join(home_dir,'src'), os.path.join(home_dir,'src','utilities'),])
env.Append(CPPPATH=[os.path.join(home_dir,'src','cosmology'),])
env.Append(LIBS=['CCfits', 'tmv', 'blas', 'cfitsio', 'gsl'])
env.Append(LIBPATH = ['/vol/software/software/tools/tmv/tmv0.72/x86_64/lib', '/usr/local/lib'])
#env.Append(CPPPATH=['-stdlib=libstdc++'])

def ReadFileList(fname):
    """
    This reads a list of whitespace separated values from the input file fname
    and stores it as a list.  We will make this part of the environment so
    other SConscripts can use it
    """
    try:
        files=open(fname).read().split()
    except:
        print 'Could not open file:',fname
	sys.exit(45)
    files = [f.strip() for f in files]
    return files

# subdirectory SConscript files can use this function
env['_ReadFileList'] = ReadFileList



# build the sub-objects
Export("env")
SConscript("src/SConscript")

# append options specific to the head node
env.Append(CPPPATH=['utilities',])

# specify the sub-objects
sub_objects = '''
	    src/cosmology/.obj/Cosmology.o
	    src/cosmology/.obj/GrowthFunction.o
	    src/cosmology/.obj/HuEisenstein.o
	    src/cosmology/.obj/LensQuantities.o
	    src/cosmology/.obj/LinearPowerSpectrum.o
	    src/cosmology/.obj/Bispectra.o
	    src/utilities/.obj/StringStuff.o
	    src/utilities/.obj/Bins.o
	    src/utilities/.obj/Mesh.o
	    src/utilities/.obj/GTable.o
	    src/utilities/.obj/odeint.o
	    src/utilities/.obj/Matrix.o
	    src/.obj/GGLens.o
	    src/.obj/LensObjects.o
	    src/.obj/StarMaskObjects.o
	    src/.obj/GAMAObjects.o
	    src/.obj/RandomObjects.o
	    src/.obj/GAMARandomObjects.o
	    src/.obj/SourceObjects.o
	    src/.obj/RCSLenSObjects.o
	    src/.obj/KiDSObjects.o
	    src/.obj/Shear.o
	    '''.split()

# build the main programs

#test = env.Program(target='GGLensTest', source=sub_objects+['src/GGLensDriver.cpp',])
#starhalo = env.Program(target='gglens_starhalo', source=sub_objects+['src/gglens_starhalo.cpp',])
gama_kids = env.Program(target='gglens_gama_kids', source=sub_objects+['src/gglens_gama_kids.cpp',])
gama_kids_angular = env.Program(target='gglens_gama_kids_angular', source=sub_objects+['src/gglens_gama_kids_angular.cpp',])
#random_kids = env.Program(target='gglens_random_kids', source=sub_objects+['src/gglens_random_kids.cpp',])
#env.Install('bin', [test, gama_kids, random_kids, starhalo])
#env.Install('bin', [gama_kids])
#env.Install('bin', [random_kids])
#env.Alias('install', 'bin')

# (eventually, build a library code)
