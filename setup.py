import osmipy

from setuptools import setup
try:  # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError:  # for pip <= 9.0.3
    from pip.req import parse_requirements

pkgs = []
dependency_links = []
for pkg in parse_requirements('requirements.txt', session=False):
    if pkg.link:
        dependency_links.append(str(pkg.link))
    else:
        pkgs.append(str(pkg.req))

setup(
    name=osmipy.__name__,
    packages=['osmipy'],
    version=osmipy.__version__,
    author=osmipy.__author__,
    author_email=osmipy.__email__,
    description=osmipy.__doc__,
    classifiers=[
        'Environment :: Scientific',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3'
    ],
    install_requires=pkgs,
    test_suite='tests'
)