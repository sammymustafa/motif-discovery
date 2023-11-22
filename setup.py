from setuptools import setup, find_packages

setup(
    name='motif-discovery',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'matplotlib',
        'seaborn',
        'pyranges',
        'biopython',
        'pysam',
        'pyfaidx',
        'logomaker',
        'anndata',
        'torch',
        'scikit-learn',
        'tqdm',
        'hmmlearn'
    ],
    python_requires='>=3.6',
)
