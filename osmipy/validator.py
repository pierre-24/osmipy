from osmipy import visitor


class ValidationException(Exception):
    pass


class Validator(visitor.ASTVisitor):
    """Visitor that:

    + Checks that ring bonds are valid:

      - They must go by pair ;
      - An atom cannot bond to itself ;
      - An atom cannot bond to a direct neighbour via a ring bond ;
      - The same ring bond cannot be defined twice ;
      - The same bond symbol must be used, if both are given.

      If everything is ok, then the ``target`` parameter of the ``Ringbond`` is set.
    + Keeps a dictionary of the atoms id, ``atom_ids``, and check that they are all uniques.

    :param node: the node
    :type node: Chain
    """
    def __init__(self, node):
        super().__init__(node)

        self.ring_ids = {}
        self.atom_ids = {}
        self.next_id = 0

    def validate(self, shift_id=0):
        """Validate the AST

        :param shift_id: shift all atom id
        :type shift_id: int
        """
        if self.node is not None:
            self.ring_ids = {}
            self.atom_ids = {}
            self.next_id = 0

            self._start(shift_id=shift_id)

            all_pairs = []

            for i, n in self.ring_ids.items():
                if len(n) % 2 != 0:
                    raise ValidationException('ring id {}: not used an even number of time'.format(i))
                for k in range(int(len(n) / 2)):
                    rb1 = n[k * 2]
                    rb2 = n[k * 2 + 1]

                    if rb1.bond is not None and rb2.bond is not None:
                        if rb1.bond != rb2.bond:
                            raise ValidationException(
                                'ring_id {}: bond are not the same ({} != {})'.format(i, n[k * 2], n[k * 2 + 1]))

                    if rb1.parent == rb2.parent:
                        raise ValidationException('ring id {}: bond to same atom'.format(i))

                    if rb1.parent.parent.right == rb2.parent.parent:
                        raise ValidationException('ring id {}: bond to a direct pair'.format(i))

                    rb1.target, rb2.target = rb2.parent, rb1.parent
                    pair = (rb1.parent, rb2.parent)
                    if id(rb2.target) < id(rb1.target):
                        pair = tuple(reversed(pair))

                    if pair in all_pairs:
                        raise ValidationException('ring id {}: this ring bond was already defined!'.format(i))
                    else:
                        all_pairs.append(pair)

    def visit_atom(self, node, shift_id=0, *args, **kwargs):
        """Just update the list of ids

        :param node: node
        :type node: qcip_tools.smiles.Atom
        :param shift_id: shift all atom id
        :type shift_id: int
        """

        super().visit_atom(node, *args, **kwargs)

        if node.atom_id > -1:
            if shift_id != 0:
                node.atom_id += shift_id

            if node.atom_id not in self.atom_ids:
                self.atom_ids[node.atom_id] = node
                if self.next_id <= node.atom_id:
                    self.next_id = node.atom_id + 1
            else:
                raise ValueError('two atoms share the same id: {}!'.format(node.atom_id))

    def visit_ringbond(self, node, *args, **kwargs):
        """All the magic happen here: count number of times when the ring_id was used, and the type of bond

        :param node: node
        :type node: RingBond
        :rtype: list of int
        """

        super().visit_ringbond(node, *args, **kwargs)

        if node.ring_id in self.ring_ids:
            self.ring_ids[node.ring_id].append(node)
        else:
            self.ring_ids[node.ring_id] = [node]

    def update(self, fragment):
        """Update the validator with a new fragment

        :param fragment: the new fragment
        :type fragment: Chain
        """

        validator = Validator(fragment)

        validator.validate(shift_id=self.next_id)
        self.atom_ids.update(validator.atom_ids)
        self.next_id = validator.next_id

        for i, j in validator.ring_ids.items():
            if i in self.ring_ids:
                self.ring_ids[i].extend(j)
            else:
                self.ring_ids[i] = j
