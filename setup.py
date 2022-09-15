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
    description="Chisel package part of the Poseydon suite",
    install_requires=['numpy>=1.20.0' , 
                      'pandas>=1.3.0', 
                      'pysam>=0.18.0', 
                      'pybigwig', 
                      'biopython>=1.79', 
                      'joblib>=1.1.0', 
                      'matplotlib>=3.5.0', 
                      'denseweight==0.1.2'],
    license="MIT license",
    include_package_data=True,
    keywords='poseydon_chisel',
    name='poseydon-chisel',
    packages=['poseydon_chisel'],
    url='https://github.com/Husam94/poseydon-chisel',
    version='0.0.3',
    zip_safe=False
)
