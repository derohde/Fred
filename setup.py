import os
import re
import sys
import sysconfig
import platform
import subprocess

from glob import glob
from distutils.version import LooseVersion
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools

this_dir = os.path.dirname(os.path.abspath(__file__))

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=this_dir):
        Extension.__init__(self, name, 
                           include_dirs=[os.path.join(this_dir, 'include'), os.path.join(this_dir, 'pybind11')],
                           sources=sorted(glob("src/*.cpp")))
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '.', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                         out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
            
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(),
                extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.call(['git', 'init'])
        subprocess.call(['git', 'submodule', 'add', 'https://github.com/pybind/pybind11.git'])
        subprocess.call(['git', 'submodule', 'update', '--init', '--recursive'])
        subprocess.check_call(['cmake', "{}".format(ext.sourcedir)] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)

setup(
    name='Fred-Frechet',
    version='1.7.4',
    author='Dennis Rohde',
    author_email='dennis.rohde@tu-dortmund.de',
    description='A fast, scalable and light-weight C++ Fréchet distance library, exposed to python and focused on (k,l)-clustering of polygonal curves.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="http://fred.dennisrohde.work",
    packages=setuptools.find_packages(),
    ext_package="Fred",
    ext_modules=[CMakeExtension('backend')],
    install_requires=['cvxopt', 'matplotlib'],
    cmdclass=dict(build_ext=CMakeBuild),
    data_files=[('pybind11', ['.gitmodules'])],
    zip_safe=False,
)
