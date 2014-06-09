#
# SConstruct
# 

# CPPPATH: sets the include path for C++
env = Environment()
env.Append(CPPPATH=['utilities',])
env.Append(LIBS=[])

# build the sub-objects
Export("env")
SConscript("utilities/SConscript")

# specify the sub-objects
sub_objects = '''
	    utilities/StringStuff.o
	    utilities/Bins.o
	    LensObjects.o
	    Mesh.o
	    '''.split()

# build the objects
env.Object('LensObjects.cpp')
env.Object('Mesh.cpp')

# build the main programs
env.Program(target='GGLensTest', source=sub_objects+['GGLensDriver.cpp',])

# (eventually, build a library code)