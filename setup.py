#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'astropy',
    'ligo-gracedb',
    'ligo.skymap',
    'numpy',
    'matplotlib',
    'healpy',
]

test_requirements = ['pytest>=3', ]

setup(
    author="VTDA Group",
    author_email='karthik_yadavalli@g.harvard.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Pipeline to process LIGO alerts for the LIGHETR collaboration",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='lighetr_alert_system',
    name='lighetr_alert_system',
    packages=find_packages(include=['lighetr_alert_system', 'lighetr_alert_system.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/sky5265/LIGHETR_Alert_System',
    version='0.1.0',
    zip_safe=False,
)
