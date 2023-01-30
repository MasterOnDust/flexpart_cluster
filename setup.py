import setuptools 
from numpy.distutils.core import Extension
from setuptools import find_packages

ext1 = Extension(name='fpcluster.trajectory_clustering', 
            sources=['fpcluster/clustering/trajectory_clustering.pyf', 'fpcluster/clustering/trajectory_clustering.f90', 
                    'fpcluster/clustering/centerofmass.f90'])

if __name__ =="__main__":
    from numpy.distutils.core import setup
    
    setup(name='fpcluster',
    version="1.0.1",
    description= 'Python wrapper of FLEXPART trajectory clustering subroutine',
    author ="Ove Haugvaldstad",
    author_email='ovehaugv@outlook.com',
    packages=find_packages(),
    ext_modules = [ext1])