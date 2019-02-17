from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy
# from Cython.Build import cythonize

try:
    # from Cython.Distutils.extension import Extension
    from Cython.Distutils import build_ext
except ImportError:
    from setuptools import Extension
    USING_CYTHON = False
else:
    USING_CYTHON = True

ext = '.pyx' if USING_CYTHON else '.c'

extensions = [
    Extension(
        "flashgamma.gamma_evaluation_c",
        ["flashgamma/gamma_evaluation_c" + ext],
        include_dirs=[numpy.get_include()],
    ),
]

if USING_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="flashgamma",
    version="0.0.1",
    author="Samuel Peet",
    author_email="samuel.peet.physics@gmail.com",
    description="Analysis of radiotherapy dose distributions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samuelpeet/flashgamma",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    ext_modules = extensions
)