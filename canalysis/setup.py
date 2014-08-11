from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

extensions = [Extension("int_scattering", ["int_scattering.pyx"]),
           Extension("particle_distance", ["particle_distance.pyx"]) ] 
#setup(
#            cmdclass = {'build_ext': build_ext},
#                ext_modules = cythonize(extensions, gdb_debug=True)
#                )
setup(
            cmdclass = {'build_ext': build_ext},
                ext_modules = extensions
                )
