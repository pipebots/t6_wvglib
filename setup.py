import pathlib
from setuptools import setup, find_packages


HERE = pathlib.Path(__file__).parent

README = (HERE / "README.md").read_text()

setup(
    name="wvglib",
    version="0.0.0",
    description="Metal and dielectric waveguide calculations",
    long_description=README,
    long_description_content_type="text/markdown",
    url="",
    author="Viktor Doychinov",
    author_email="eenvdo@leeds.ac.uk",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
    packages=find_packages(exclude=("tests",)),
    include_package_data=True,
    install_requires=["numpy", "scipy", "rflib"],
)
