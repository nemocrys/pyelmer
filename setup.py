# https://packaging.python.org/tutorials/packaging-projects/
import versioneer
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyelmer",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Arved Enders-Seidlitz",
    author_email="arved.enders-seidlitz@ikz-berlin.de",
    description="A python interface to Elmer.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nemocrys/pyelmer",
    packages=["pyelmer", "pyelmer.test"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "gmsh",
        "pyyaml",
        "matplotlib",
        "numpy",
        "objectgmsh",
        "pandas",
    ],
    python_requires=">=3.7",
)
