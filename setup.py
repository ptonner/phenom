from setuptools import setup, find_packages
from os import path
from codecs import open

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='phenom',

    version='0.0.1a',

    description='A hierarchical non-parametric microbial phenotype model',
    long_description=long_description,

    url='https://github.com/ptonner/phenom',

    author='Peter Tonner',
    author_email='peter.tonner@gmail.com',

    license='Apache 2.0',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',

                'Environment :: Console',

        'Intended Audience :: Developers',
                'Intended Audience :: Science/Research',
        'Topic :: Database :: Database Engines/Servers',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ],

    keywords='statistics microbiology',

    packages=find_packages(exclude=['data', 'tests']),
    include_package_data=True,
    package_data={'phenom': ['stan/*.stan']},

    install_requires=['pandas', 'matplotlib',
                      'numpy', 'pystan', 'patsy',
                      'scipy', 'attrs'],
    test_requires=["pystest", "pytest-runner"],

)
