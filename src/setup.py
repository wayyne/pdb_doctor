import os
from setuptools import setup, find_packages

# Read the long description from the README file if it exists.
long_description = ""
readme_path = os.path.join(os.path.dirname(__file__), "README.md")
if os.path.exists(readme_path):
    with open(readme_path, "r", encoding="utf-8") as fh:
        long_description = fh.read()

setup(
    name="pdbdr",
    version="1.0.0",
    description=(
        "A toolkit for processing, aligning, and stitching protein "
        "structures."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Your Name",  # Replace with your name
    author_email="your.email@example.com",  # Replace with your email
    url="https://github.com/yourusername/pdbdr",  # Optional: Update URL if available
    packages=find_packages(),
    install_requires=[
        "numpy",
        "torch",
        "tqdm",
        "requests",
        "esm",
    ],
    python_requires=">=3.7",
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
)

