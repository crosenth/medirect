import setuptools
import sys

setuptools.setup(
    author='Chris Rosenthal',
    author_email='crosenth@gmail.com',
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Operating System :: OS Independent',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only'],
    description='multiprocessed ncbi edirect and ftract',
    entry_points={
        'console_scripts': [
            'mefetch=medirect.mefetch:run',
            'ftract= medirect.ftract:run']},
    install_requires=['biopython>=1.68', 'retrying>=1.3.3'],
    keywords=['ncbi', 'edirect', 'multiprocessing',
              'entrez', 'bioinformatics'],
    license='GPLv3',
    name='medirect',
    packages=setuptools.find_packages(exclude=['tests']),
    python_requires='>=3.4',
    version='0.21.0',
    url='https://github.com/crosenth/medirect',
)
