import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.verbose = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


def extract_version(module='utide'):
    version = None
    fdir = os.path.dirname(__file__)
    fnme = os.path.join(fdir, module, '__init__.py')
    with open(fnme) as fd:
        for line in fd:
            if (line.startswith('__version__')):
                _, version = line.split('=')
                # Remove quotation characters.
                version = version.strip()[1:-1]
                break
    return version

rootpath = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    return open(os.path.join(rootpath, *parts), 'r').read()

long_description = '{}\n{}'.format(read('README.rst'), read('LICENSE.txt'))

with open('requirements.txt') as f:
    require = f.readlines()
install_requires = [r.strip() for r in require]

setup(name='UTide',
      version=extract_version(),
      license='MIT',
      long_description=long_description,
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Education',
                   'License :: OSI Approved :: MIT License',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Education',
                   ],
      description='Python distribution of the MatLab package UTide.',
      url='https://github.com/wesleybowman/UTide',
      platforms='any',
      keywords=['oceanography', 'tides'],
      install_requires=install_requires,
      packages=['utide'],
      package_data={'utide': ['data/*.npz']},
      tests_require=['pytest'],
      cmdclass=dict(test=PyTest),
      author=['Wesley Bowman'],
      author_email='wesley.bowman23@gmail.com',
      zip_safe=False)
