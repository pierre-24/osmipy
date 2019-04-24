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

