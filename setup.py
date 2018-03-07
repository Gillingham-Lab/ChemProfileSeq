from setuptools import setup, find_packages

build_exe_options = {
    "includes": [
        "sys",
        "subprocess",
        "os",
        "numpy",
        "matplotlib",
        "blessings",
        "pysam"
    ],
    "excludes": [],
}

setup(
    name="chemProfileSeq",
    version="1.0",
    description="Toolchain to analyze chemical profiles in RNA.",
    url="https://github.com/Gillingham-Lab/ChemProfileSeq",

    author="Basilius Sauter",
    author_email="basilius.sauter@unibas.ch",

    packages=find_packages(),
    install_requires=[
        "argh",
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
        "blessings",
        "pysam"
    ],

    entry_points={
        "console_scripts": [
            'chemProfileSeq=chemProfileLib.cli.main_app:main',
        ],
    },
    #scripts=["bin/chemProfileSeq"]
)