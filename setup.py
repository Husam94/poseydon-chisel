from setuptools import setup

setup(
    author = "Husam Abdulnabi",
    author_email = 'husam.abdulnabi@gmail.com',
    python_requires = '>=3.7',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education', 
        'Intended Audience :: Developers', 
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering' 
    ],
    description="Chisel package part of the Poseidon suite",
    install_requires=['numpy' , 'pandas', 'random', 'math', 'pickle', 'pysam', 
                      'pyBigWig', 'itertools', 'biopython', 'joblib', 'matplotlib', 
                      'denseweight']
    license="MIT license",
    include_package_data=True,
    keywords='poseidon_chisel',
    name='poseidon-chisel',
    packages=['poseidon-chisel']
    url='https://github.com/Husam94/poseidon-chisel',
    version='0.0.1',
    zip_safe=False
)
