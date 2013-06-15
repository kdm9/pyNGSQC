from setuptools import setup

desc = """
pyNGSQC is a python re-write of many of the functions of the fastx tookit,
seqtk, FastQC and other similar next-gen sequencing quality control tools

It is designed to be fast (and particularly, compatible with PyPy), memory
efficient, and import-able into existing or new python programs or pipelines.
"""

setup(
    name="pyNGSQC",
    packages=['pyngsqc', ],
    version="0.1a",
    description=desc,
    author="Kevin Murray",
    author_email="k.d.murray.91@gmail.com",
    url="https://github.com/kdmurray91/pyNGSQC",
    keywords=["http", "multipart", "post", "urllib2"],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later " +
            "(GPLv3+)",
        ],
    test_suite="test",
    )
