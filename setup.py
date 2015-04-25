import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand

rootpath = os.path.abspath(os.path.dirname(__file__))


class PyTest(TestCommand):
    """python setup.py test"""
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--strict', '--verbose', '--tb=long', 'tests']
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


def readme():
    with open('README.md') as f:
        return f.read()


def extract_version():
    version = None
    fname = os.path.join(rootpath, 'utide', '__init__.py')
    with open(fname) as f:
        for line in f:
            if (line.startswith('__version__')):
                _, version = line.split('=')
                version = version.strip()[1:-1]  # Remove quotation characters
                break
    return version

setup(name='UTide',
      version=extract_version(),
      description='Python distribution of the MatLab package UTide.',
      long_description=readme(),
      url='https://github.com/wesleybowman/UTide',
      author='Wesley Bowman',
      author_email='wesley.bowman23@gmail.com',
      maintainer='Wesley Bowman',
      license='MIT',
      packages=['utide', 'utide/tests'],
      package_data={'utide': ['data/*.mat']},
      tests_require=['pytest'],
      extras_require=dict(testing=['pytest']),
      cmdclass=dict(test=PyTest),
      zip_safe=False)
