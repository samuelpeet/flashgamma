import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="flashgamma",
    version="0.0.2",
    author="Samuel Peet",
    author_email="samuel.peet.physics@example.com",
    description="Analysis of radiotherapy dose distributions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samuelpeet/flashgamma",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'scipy'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)