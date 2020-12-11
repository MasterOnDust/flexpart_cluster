from numpy.distutils.core import Extension
from setuptools import find_packages

ext1 = Extension(name='fpcluster.clustering', 
            sources=['src/fpcluster/fpcluster.pyf', 'src/fpcluster/clustering.f90','src/fpcluster/par_mod.f90'])
ext2 = Extension(name='fpcluster.clustering_whole',
            sources=['src/fpcluster/fpcluster_whole.pyf','src/fpcluster/clustering_whole.f90','src/fpcluster/par_mod.f90'])

if __name__ =="__main__":
    from numpy.distutils.core import setup
    
    setup(name='fpcluster',
    version="0.0.1",
    description= 'Python wrapper of FLEXPART trajectory clustering subroutine',
    author ="Ove Haugvaldstad",
    author_email='ovehaugv@outlook.com',
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    ext_modules = [ext1,ext2])