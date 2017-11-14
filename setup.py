from distutils.core import setup, Extension

boost_module = Extension('_seqbdd',
                         sources=['src/seqbdd/_seqbdd.cpp',
                                  'src/seqbdd/matrix/matrix.cpp'],
                         include_dirs=['/usr/local/include'],
                         library_dirs=['/usr/local/lib'],
                         libraries=['boost_python3'],
                         extra_compile_args=['-std=c++14', '-O2', '-Wall'])

setup(name='seqbdd',
      version='0.1.0',
      author='daiwahome',
      package_dir={'': 'src'},
      packages=['seqbdd'],
      ext_modules=[boost_module])
