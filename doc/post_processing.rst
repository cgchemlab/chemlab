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

It allows changing property at certain topological distance from the reacted particle..

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


Remove neighbour bond
----------------------------

Remove a bond that is at certain topological distance from the reacted particle.

.. glossary::

   ext_type
      RemoveNeighboursBonds

   bonds_to_remove
      The comma separated list of ``bonds_to_remove`` definitions.
      The syntax of single ``bonds_to_remove``:
      ``anchor_type -> type1:type2:distance``
      where ``anchor_type`` is the type of *root* particle, from where the distance is measured,
      ``type1:type2`` is the bond to remove (in terms of particle types) at ``distance``.


Release molecule
----------------------------

Remove the distance constraint and release a molecule linked to host molecule.

.. glossary::

   ext_type
      ReleaseMolecule

   host_type
      The type of host particle.

   target_type
      The type of dummy non-released particle.

   eq_length
      The distance at which the dummy particle is placed from the host particle.

   alpha
      The rate constant used to fade in the dummy particle after released.

   init_res
      The initial resolution of the dummy particle after released.

   final_type
      The final type of the dummy particle after completely fade in (resolution 1)

   cache_file
      The cache file for dummy particles (optional)

   replicate
      How many dummy particles should be attached to a single host particle (default: 1)

   release_on

      - ``type`` then the particle will be released whenever the host particle change a type (default)
      - ``release_on`` remove whenever host particle react with other particle

   release_count
      number of particles to release (default: 1)


Join molecule
----------------------------

Add the distance constraint and make invisible one of the particles.

.. glossary::

   ext_type
      JoinMolecule


Freeze regions
----------------------------

Define that the box edges will play a role a freeze regions. Whenever a particle of given
type reaches this region, it becomes invisible.

.. glossary::

   ext_type
      FreezeRegion


Change randomly particle type
----------------------------

Chang randomly particle types during the simulation

.. glossary::

   ext_type
      ChangeParticleType


Simulate ATRP class of the reactions
----------------------------

Chang randomly particle types during the simulation

.. glossary::

   ext_type
      ATRPActivator