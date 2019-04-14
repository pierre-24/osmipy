Installing ``osmipy``
=====================

Normal install
--------------

To just use ``qcip_tools`` in your Python projects, simply use pip:

.. code-block:: bash

    pip3 install git+ssh://git@github.com:pierre-24/osmipy.git@dev

Note that ``--user`` can allow you to install the package without being superuser (see `here <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_).
On the other hand, you can install it in a *virtualenv*.

You can also add it to your ``requirements.txt`` or Pipenv's  ``Pipfile``:

.. code-block:: text

    # requirements.txt style:
    git+ssh://git@github.com:pierre-24/osmipy.git@dev

    # Pipfile style
    qcip-tools = {ref = "dev", git = "ssh://git@github.com:pierre-24/osmipy.git"}


Installation for contributors
-----------------------------

To contribute to the project, you need to clone the repository:

+ Clone it: ``git clone git@github.com:pierre-24/osmipy.git``.
+ Install pipenv: ``pip3 install pipenv``
+ Install virtualenv and dependencies: ``make init``.

You can launch the tests series with ``make test``

Don't forget to check the `contribution rules <contributing.html>`_.