import os
from .cov import covariance, correlation

with open(f"{os.path.dirname(os.path.abspath(__file__))}/VERSION", "r") as f:
    __version__ = f.read().strip()
