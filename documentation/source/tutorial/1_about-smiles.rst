Parsing SMILES
==============

About SMILES
------------

SMILES (simplified molecular input entry system) is a "chemical notation system designed for modern chemical information processing" invented by D. Weininger and then developed by Daylight Chemical Information Systems (simply referred as "Daylight" is most websites).
There is also `an open standard <http://opensmiles.org/opensmiles.html>`_ specification that will be followed here.

A SMILES string is basically a linear representation of a molecule created from a traversal of the molecular graph (in its simplest version, a labeled graph whose vertices are the atoms and edges are the chemical bonds).
This is actually equivalent to the construction of a `spanning tree <https://en.wikipedia.org/wiki/Spanning_tree>`_.
Therefore, many SMILES are possible (including very peculiar ones) for a given molecule (depending, basically on the starting vertex, and on the edges taken during the visit).


Such SMILES string consists of a series of characters, without space.

.. note::

    Complete explanations of the SMILES syntax are provided for example on the `SMILES theory page of Daylight <http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>`_.
    A small introduction is given in the next sections.



Atoms
_____

The atomic symbol is used to specify an atom, **in square bracket**, except for the atom in the *organic subset*: B, C, N, O, S, P, Cl, Br (plus ``as`` and ``se``, but only in their aromatic form).

+ For the atom outside this organic subset, the number of hydrogen and charge must be specified, e.g. ``[H+]``, ``[Al+3]``, ``[BH3]``, ``[NaH]`` ...
+ For atom in the organic susbet, the number of hydrogen is assumed, in order to match the lowest valence, consistant with explicit bond. Thus,

  + ``C`` represent the methane ;
  + ``Cl`` represent the hydrogen chloride ;
  + ``O`` represent the water molecule ;
  + ...

  To give an explicit count of hydrogen or a charge, square bracket can be used, for example ``[NH4+]`` for the ammonium cation.
+ To differentiate atom that are in an aromatic cycle, those atoms must be lowercase.
  For example, ``c`` is used for the carbon in benzene, instead of ``C``.
  This is only valid for atom in the organic subset, except Br and Cl.


Bonds
_____

Adjacent atoms are assumed to be connected to each other by a single (or aromatic) bond.
Thus ``CCO`` is used to represent ethanol, while ``COC`` stands for dimethyl ether.
Double and triple bonds are explicitly marked, with ``=`` or ``#``, respectively.
The symbol ``-`` can be used for single bond, but is usually implicit.

For example,

+ ``CCC`` represent the propane ;
+ ``C=CC=C`` represent the butene ;
+ ``O=C=O`` stands for carbon dioxyde ;
+ ``CC#N`` represent the acetonitrile.

Branches and cycles
___________________

Branches are given in parentheses, which are linked to the previous atom (to the left).
Parentheses can be nested and stacked.

.. figure:: ../images/tuto_branches.png
    :align: center

    Different branched molecules and their corresponding SMILES.


About canonical SMILES
______________________

To handle this, people propose canonicalization schemes, that consists in two parts:

1. Renumber the atoms (the vertices of the graph) based on some invariant(s), in a way that is (normally) unique ;
2. Starts by the atom with the lowest number out of the previous step, and perform some branching decision during the exploration (once again, in a *unique* way).

The critical part is actually the first one.
For example, the algorithm developed by D. Weininger in 1989 (CANGEN) fails for some structures.

About "stereochemistry"
-----------------------

**TL;DR:**  SMILES is not about absolute stereochemistry (in the CIP sense), it is local stereochemistry.

Background
__________

Implementation of the stereo perception is directly inspired by `OpenBabel implementation of stereoconcepts <http://openbabel.org/dev-api/classOpenBabel_1_1OBStereoBase.shtml>`_, which was implemented in the same spirit as SMILES.
In SMILES (and therefore here), one only cares about **local** stereochemistry: even though a carbon may not be asymetric (in the CIP sense), it may present a **local** configuration.
There is therefore **no** correspondence between the configuration in a SMILE string and the absolute configuration (although if the vertex are given in the correct order in the CIP way, that corresponds to R/S stereo configuration for tetrahedral carbons).
One only cares about the correspondance between the vertices, given in a certain order, and the "reality" of the structure.

In practice
___________

Stereochemistry is defined for any set of 4 atoms that contains at most 1 hydrogen (non-planar) or 2 (planar).
The implementation should recognize at least two kind of stereochemistry:

.. figure:: ../images/stereo.png
    :align: center


+ Planar stereo config (square planar) ;
+ Non planar stereo config (clockwise/anticlockwise for "asymetric" carbons, allene configuration).

The purpose of this stereo implementation is to match the set of indices *a, b, c and d* (which may be given in any order) with the set *0, 1, 2 and 3*.
The value of the stereo object is the order in which the *a, b, c and d* indices must be read to match *0, 1, 2 and 3*.

For the non planar stereo configuration, there is two possible values: clockwise (``@@``) and counterclockwise (``@``).
For the planar stereo configuration, there is three possible values: ``U``, ``Z`` and ``4``.

.. figure:: ../images/stereo_vals.png
    :align: center

    The question is always "how to read ``(a,b,c,d)`` so that the sequence matches ``(0,1,2,3)``".
    To do so, the first part is to set a ``start`` atom (the others are the ``refs`` in the implementation) and set the corresponding ``value``.


Parsing SMILES with ``osmipy``
------------------------------

.. code-block:: python

    from osmipy import parse
    smiles = parse('C1(Cl)CC1')

and that's it !