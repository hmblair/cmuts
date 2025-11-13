from pathlib import Path
from setuptools import setup, find_packages

DIR = Path(__file__).parent
NAME = DIR.name
VERSION = (DIR / "VERSION").read_text().strip()
AUTHOR = "Hamish M. Blair"
EMAIL = "hmblair@stanford.edu"
URL = "https://github.com/hmblair/cmuts"
LICENSE = "MIT"

setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    url=URL,
    license=LICENSE,
    package_dir={"": "src/python"},
    packages=find_packages(where="src/python"),
)
