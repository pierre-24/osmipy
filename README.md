# osmipy (Open SMILES Python)

A simple python library to interpret SMILES ([Simplified Molecular InputLine Entry System](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)).

By [Pierre Beaujean](https://pierrebeaujean.net/)

## What ?

Currently the library is able to

+ Parse a correct SMILES string and extract an AST ;
+ Check a few extra stuffs (if ring bonds are correct and if `atoms_id` are uniques). 

Minimalistic example:

```python
from osmipy import parse
obj = parse('C12(CCCCC1)CCCCC2')

# AST starts in `obj.node`
# TODO: more stuffs
```

Check  [here](https://pierre-24.github.io/osmipy/) for more details.

## Install and/or contribute

The installation procedure is detailed [here](https://pierre-24.github.io/osmipy/install.html), while the contribution rules are listed [here](https://pierre-24.github.io/osmipy/contributing.html).
