from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='UTide',
      version='1.0_exp',
      description='Python distribution of the MatLab package UTide.',
      long_description=readme(),
      url='https://github.com/wesleybowman/UTide',
      author='Wesley Bowman',
      author_email='wesley.bowman23@gmail.com',
      license='MIT',
      packages=['utide'],
      zip_safe=False)
