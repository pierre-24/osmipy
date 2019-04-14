class NodeVisitor(object):
    """Implementation of the visitor pattern. Expect ``visit_[type](node)`` functions, where ``[type]`` is the
    type of the node, **lowercase**.
    """

    def visit(self, node, *args, **kwargs):
        method_name = 'visit_' + type(node).__name__.lower()
        visitor = getattr(self, method_name, self.generic_visit)
        return visitor(node, *args, **kwargs)

    def generic_visit(self, node):
        raise Exception('No visit_{} method'.format(type(node).__name__.lower()))


class ASTVisitor(NodeVisitor):
    """Generic visitor for the AST

    :param node: node
    :type node: Chain
    """
    def __init__(self, node):
        self.node = node

    def _start(self, *args, **kwargs):
        """Start the visit
        """
        self.visit(self.node, *args, **kwargs)

    def visit_chain(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: Chain
        """
        self.visit(node.left, *args, **kwargs)
        if node.bond is not None:
            self.visit(node.bond, *args, **kwargs)
        if node.right is not None:
            self.visit(node.right, *args, **kwargs)

    def visit_branchedatom(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: BranchedAtom
        """

        self.visit(node.atom, *args, **kwargs)

        for ringbond in node.ring_bonds:
            self.visit(ringbond, *args, **kwargs)
        for branch in node.branches:
            self.visit(branch, *args, **kwargs)

    def visit_atom(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: qcip_tools.smiles.Atom
        """
        pass

    def visit_bond(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: qcip_tools.smiles.Bond
        """
        pass

    def visit_branch(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: Branch
        """
        if node.bond is not None:
            self.visit(node.bond, *args, **kwargs)

        self.visit(node.chain, *args, **kwargs)

    def visit_ringbond(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: RingBond
        """

        if node.bond is not None:
            self.visit(node.bond, *args, **kwargs)
