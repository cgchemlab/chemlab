Post-processing modules
========================

The list of post-processing modules available directly from `.cfg` file.
Those modules have to be attached to reaction groups.

In the example below, we defined a reaction group ``abc`` with three extensions: ``release_molecule``,
``join_molecule`` and ``remove_bond``.
The name of those extensions have to corresponds to the ``ext_`` sections in the ``.cfg`` file.
Here, we only show section ``ext_remove_bond`` that corresponds to extensions ``remove_bond``.

Every *extension* section have at least two variables:

 - ``ext_type`` defines the extension type (required)
 - ``invoke_on`` declares on which reactant type the extension will be invoked.
   The possible options are: ``type_1``, ``type_2``, ``both`` (default: ``both``)

Every of the *extension* section have to set ``ext_type`` variable. Other variables in
this file are related to the particular extension.


.. code-block:: ini

   [ext_remove_bond]
   ext_type:RemoveNeighboursBonds
   bonds_to_remove=C->C:E:1
   invoke_on=type_1

   [group_abc]
   potential:Harmonic
   potential_options:K=13622.3,r0=0.256395
   extensions:release_molecule,join_molecule,remove_bond

Change neighbour property
----------------------------

It allows changing property at certain topological distance from the reactant.

.. glossary::

   ext_type
      ChangeNeighboursProperty

   type_transfers
      The comma separated list of ``type_transfer`` definitions.
      The syntax of single ``type_transfer``:
      ``old_type:distance -> new_type(opt1,opt2)``
      where ``opt1`` are optional new parameters of the particle like:

       - mass: new mass
       - q: new charge
       - state: new chemical state

      The parameters set here will override the one defined in the topology file.