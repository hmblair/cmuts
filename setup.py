from setuptools import setup, find_packages

setup(
    name="cmuts",
    version="0.1.0",
    author="Hamish M. Blair",
    author_email="hmblair@stanford.edu",
    url="https://github.com/hmblair/cmuts",
    package_dir={"": "src/python"},
    packages=find_packages(where="src/python"),
)
