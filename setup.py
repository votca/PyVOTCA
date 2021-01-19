#!/usr/bin/env python
import os

from setuptools import setup


HERE = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(HERE, 'votcapytools', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='votcapytools',
    version=version,
    description="Generic tools to interact with Votca using Python",
    long_description=readme + '\n\n',
    author="Bjoern Baumeier",
    author_email='',
    url='https://github.com/https://github.com/votca/votcapytools/votcapytools',
    packages=[
        'votcapytools',
    ],
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='votcapytools',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
    ],
    install_requires=["numpy", "h5py"],
    extras_require={
        'test': ['coverage', 'mypy', 'pycodestyle', 'pytest>=3.9',
                 'pytest-asyncio', 'pytest-cov', 'pytest-mock'],
        'docs': ['sphinx', 'sphinx_rtd_theme']

)
