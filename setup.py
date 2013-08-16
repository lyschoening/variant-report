# encoding: utf-8
from distutils.core import setup

with open('requirements.txt', 'r') as f:
    requirements = f.read().split("\n")

setup(
    name='variant-report',
    version='0.2.0',
    packages=['annotation', 'annotation.utils'],
    url='',
    license='',
    author=u'Lars Schöning',
    entry_points={
        'console_scripts': [
            'variant-report = annotation.report:main',
            'coverage-report = annotation.coverage_report:main',
            'interval-report = annotation.interval_coverage_report:main'
        ]
    },
    author_email='lars@lyschoening.de',
    description='',
    install_requires=requirements
)
