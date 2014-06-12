#
# SConstruct
# 

# CPPPATH: sets the include path for C++
env = Environment()
env.Append(CPPPATH=['utilities','/vol/software/software/tools/tmv/tmv0.72/x86_64/include'])
env.Append(LIBS=['tmv', 'blas',])
env.Append(LIBPATH = ['/vol/software/software/tools/tmv/tmv0.72/x86_64/lib', '/usr/local/lib'])

# build the sub-objects
Export("env")
SConscript("utilities/SConscript")

# specify the sub-objects
sub_objects = '''
	    utilities/StringStuff.o
	    utilities/Bins.o
	    utilities/Mesh.o
	    GGLens.o
	    LensObjects.o
	    SourceObjects.o
	    Shear.o
	    '''.split()

# build the objects
env.Object('GGLens.cpp')
env.Object('LensObjects.cpp')
env.Object('SourceObjects.cpp')
env.Object('Shear.cpp')

# build the main programs
test = env.Program(target='GGLensTest', source=sub_objects+['GGLensDriver.cpp',])
env.Install('bin', [test,])
env.Alias('install', 'bin')

# (eventually, build a library code)
