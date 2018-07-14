from setuptools import setup, find_packages
import sys

# require
if sys.version_info < (3, 4):
    raise Exception('demo requires Python 3.4 or higher.')

def readme():
    with open('README.md', 'r') as f:
        return f.read()

def license():
    with open('LICENSE', 'r') as f:
        return f.read()

setup(
    name='goldclip',
    version='0.0.1',
    description='Analysis pipeline for GoldCLIP datasets',
    long_description=readme(),
    author='Ming Wang',
    author_email='wangm08@hotmail.com',
    keywords='GoldCLIP',
    url='http://github.com/bakerwm/goldclip',
    license=license(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: Linux',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Bash',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities',
    ],
    packages=find_packages(), # 
    install_requires=[
        'markdown',
        ],
    # data_files=[('bin', ['bin/*.sh'])],
    # scripts=['bin/funniest-joke'],
    entry_points={
        'console_scripts': ['goldclip=goldclip.main:main'],
    },
    test_suite='nose.collector',
    tests_require=['nose'],
    include_package_data=True,
    zip_safe=False
    )
