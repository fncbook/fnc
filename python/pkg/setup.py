from setuptools import setup, find_packages

setup(
    name='fncbook',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
    ],
    author='Tobin Driscoll',
    author_email='driscoll@udel.edu',
    description='Python codes for the book "Fundamentals of Numerical Computation"',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/fncbook/fncbook.py',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)