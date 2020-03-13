from pathlib import PosixPath as Path
from setuptools import setup, find_packages, Extension

my_name = Path(__file__).resolve()

here = lambda f: str(my_name.with_name(f))

with open(here("README.md"), "r") as fh:
    long_description = fh.read()

setup(
    name='AuthentiCT',
    version='0.1.0',
    packages=find_packages(),
    url='https://github.com/StephanePeyregne/AuthentiCT',
    license='GPLv3.0',
    author='Stephane Peyregne',
    author_email='stephane_peyregne@eva.mpg.de',
    description='A method to estimate the proportion of present-day DNA contamination among ancient sequences generated from single-stranded libraries',
    long_description=long_description,
    long_description_content_type="text/markdown",
    setup_requires=['Cython'],
    install_requires=['scipy', 'numpy', 'numdifftools', 'pandas'],
    python_requires='>=3.6',
    entry_points={'console_scripts': ['AuthentiCT=AuthentiCT.cli:main']},
    ext_modules=[Extension("AuthentiCT.contamination_est", ["AuthentiCT/contamination_est.pyx"], language="c++")]
)
