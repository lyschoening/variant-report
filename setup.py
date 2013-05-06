# encoding: utf-8
from distutils.core import setup

with open('requirements.txt', 'r') as f:
    requirements = f.read().split("\n")

setup(
    name='variant-report',
    version='0.0.3',
    packages=['annotation'],
    url='',
    license='',
    author=u'Lars Schöning',
    entry_points={
        'console_scripts': [
            'variant-report = annotation.report:main'
        ]
    },
    author_email='lars@lyschoening.de',
    description='',
    install_requires=requirements
)
