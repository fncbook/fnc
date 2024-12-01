from setuptools import setup, find_packages

setup(
    name='FNC',
    version='0.2.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
    ],
    author='Tobin Driscoll',
    author_email='driscoll@udel.edu',
    description='Python codes for the book "FUNDAMENTALS OF NUMERICAL COMPUTATION"',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/fncbook/fundamentals-numerical-computation',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)