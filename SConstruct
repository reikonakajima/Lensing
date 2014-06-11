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
	    utilities/Mesh.o
	    LensObjects.o
	    '''.split()

# build the objects
env.Object('LensObjects.cpp')

# build the main programs
test = env.Program(target='GGLensTest', source=sub_objects+['GGLensDriver.cpp',])
env.Install('bin', [test,])
env.Alias('install', 'bin')

# (eventually, build a library code)