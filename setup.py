"""
Biostar Recipes
"""
from setuptools import setup

from recipes import VERSION

setup(
    name="Biostar Recipes",
    version=VERSION,
    packages=[
        "recipes",
    ],
)
