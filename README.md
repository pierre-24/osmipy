# osmipy (Open SMILES Python)

A simple python library to interpret SMILES ([Simplified Molecular InputLine Entry System](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)).

By [Pierre Beaujean](https://pierrebeaujean.net/)

## What ?

Currently the library is able to

+ Parse a correct SMILES string and extract an AST ;
+ Check a few extra stuffs (if ringbonds are correct and if `atoms_id` are uniques). 

Minimal example:

```python
from osmipy import smiles_parser
obj = smiles_parser.parse('C12(CCCCC1)CCCCC2')

# AST starts in `obj.node`
# TODO: more stuffs
```

## Install and/or contribute

The installation procedure is detailed [here](https://pierre-24.github.io/osmipy/install.html), while the contribution rules are listed [here](https://pierre-24.github.io/osmipy/contributing.html).
