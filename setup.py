import setuptools
import sys

if sys.version_info < (3, 0):
    raise EnvironmentError('Please install using pip3 or python3')

setuptools.setup(
    author='Chris Rosenthal',
    author_email='crosenth@gmail.com',
    description='multiprocessed ncbi edirect and ftract',
    name='medirect',
    keywords=['ncbi', 'edirect', 'multiprocessing',
              'entrez', 'bioinformatics'],
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'mefetch=medirect.mefetch:run',
            'ftract= medirect.ftract:run']},
    version='0.5.0',
    url='https://github.com/crosenth/medirect',
    license='GPLv3',
    install_requires=['biopython>=1.68', 'retrying>=1.3.3'],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Operating System :: OS Independent',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only'])
