import os
from setuptools import setup

import versioneer


rootpath = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    return open(os.path.join(rootpath, *parts), 'r').read()


with open('requirements.txt') as f:
    require = f.readlines()
install_requires = [r.strip() for r in require]

setup(name='UTide',
      version=versioneer.get_version(),
      license='MIT',
      long_description=read('README.rst'),
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
      cmdclass=versioneer.get_cmdclass(),
      author=['Wesley Bowman'],
      author_email='wesley.bowman23@gmail.com',
      zip_safe=False)
