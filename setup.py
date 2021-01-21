#!/usr/bin/env python
import os

from setuptools import find_packages, setup

HERE = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(HERE, 'PyVOTCA', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='PyVOTCA',
    version=version,
    description="Generic tools to interact with Votca using Python",
    long_description=readme + '\n\n',
    author="Bjoern Baumeier",
    author_email='',
    url='https://github.com/votca/PyVOTCA/PyVOTCA',
    packages=find_packages(),
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
    entry_points={
        'console_scripts': [
            'xtp_gradient=PyVOTCA.xtp_gradient:main',
        ]
    },

    install_requires=["h5py", "matplotlib", "numpy", "pyyaml"],
    extras_require={
        'test': ['coverage', 'mypy', 'pycodestyle', 'pytest>=3.9',
                 'pytest-asyncio', 'pytest-cov', 'pytest-mock'],
        'docs': ['sphinx', 'sphinx_rtd_theme']
    }
)
