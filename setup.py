import os
from setuptools import setup

rootpath = os.path.abspath(os.path.dirname(__file__))


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
      packages=['utide'],
      package_data={'utide': ['data/*.mat']},
      zip_safe=False)
