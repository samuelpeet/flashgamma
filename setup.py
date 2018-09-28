from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

with open("README.md", "r") as fh:
    long_description = fh.read()

extensions = [
    Extension(
        "flashgamma.gamma_evaluation_c",
        ["flashgamma/gamma_evaluation_c.pyx"],
        include_dirs=[numpy.get_include()],
    ),
]

setup(
    name="flashgamma",
    version="0.0.2",
    author="Samuel Peet",
    author_email="samuel.peet.physics@example.com",
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
    ext_modules = cythonize(extensions, gdb_debug=True)
)