"""
A simple python library to interpret SMILES
"""

__name__ = 'osmipy'
__version__ = '0.1a'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

from osmipy.smiles_parser import parse


def tree_view(node):
    from osmipy.smiles_visitors import TreeView
    return TreeView(node)()
