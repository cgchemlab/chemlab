Functionals
===========

Nonbonded potentials
--------------------

Lennard-Jones
+++++++++++++

.. _lj:

.. math::

   U(r_{ij}) = 4\epsilon\left [ \left ( \frac{\sigma}{r_{ij}} \right)^{12} - \left ( \frac{\sigma}{r_{ij}} \right)^6 \right ]


Tabulated (conversion)
++++++++++++++++++++++

.. _tc:

The appropriate table is selected, based on current chemical conversion.


Tabulated (mixed, conversion)
+++++++++++++++++++++++++++++

Arithmetic mixing of two tabulated potentials

**Input**

 - tab1: first tabulated potential
 - tab2: second tabulated potential
 - type: the type of particles to count in order to calculate conversion
 - total_number: the total expected number of particles of given type **M**

The conversion is defined as:

.. math::

   \Phi = \frac{N_{type}}{M}

And the effective potential:

.. math::

   U(r_{ij}) = \Phi U^{tab1} + (1-\Phi) U^{tab2}


Bonded potentials
-----------------

Harmonic bond
+++++++++++++

.. _eq1:

.. math::

   U(r) = K(r-r0)^2

FENE bond
++++++++++++++++

.. _eqFENE:

.. math::

   U(r) = -\frac{1}{2} K b^2 log \left( 1 - \frac{r^2}{b^2} \right)

Harmonic angle
++++++++++++++

.. _eq2:

.. math::

   U(\theta) = K(\theta - \theta_0)^2


Cosine angle
++++++++++++

.. _eq3:

.. math::

   U(\theta) = K(1.0 + cos(\theta - \theta_0))

Harmonic n-cosine dihedral
++++++++++++++++++++++++++

.. _eq4:

.. math::

   U(\phi) = K(1 + cos(multiplicity*\phi - \phi_0));


Ryckaert Bellemans dihedral
+++++++++++++++++++++++++++

.. _eq5:

.. math::

   U(\phi) = \sum^{5}_{n=0} K_n cos^n(\phi)


Dihedral Harmonic
++++++++++++++++++++++++++++

.. _eq7:

.. math::
   
   U(\phi) = \frac{1}{2} K (\phi - \phi_0)^2




Topology file
-------------

[ bondtypes ]
+++++++++++++

========================  =====  =======
Name of interaction       func   params
========================  =====  =======
Harmonic eq1_             1      r0, K [1]_
FENE eqFENE_              7      r0, K
Tabulated                 8      table index
========================  =====  =======

.. [1] Force constant internally divided by 2.0

[ angletypes ]
++++++++++++++

========================  =====  =======
Name of interaction       func   params
========================  =====  =======
Harmonic eq2_             1      theta0 (deg), K [2]_
Tabulated                 8      table index
Cosine   eq3_             11     theta0 (deg), K [2]_
========================  =====  =======

.. [2] Force constant internally divided by 2.0

[ dihedraltypes ]
+++++++++++++++++

========================  =====  =======
Name of interaction       func   params
========================  =====  =======
HarmonicNCos  eq4_        1      phi0 (deg), K, multiplicity
Ryckaert Bellemans  eq5_  3      K0, K1, K2, K3, K4, K5
Tabulated                 8      table index
Harmonic  eq7_            12     phi0 (deg), K
========================  =====  =======


[ nonbond_params ]
++++++++++++++++++

==============================  ====  ======
Name of interaction             func  params
==============================  ====  ======
Lennard-Jones       lj_         1     sigma*, epsilon*
Tabulated                       8     filename*
Tabulated (conversion) tc_      9     filename*, type, total number, p_min, p_max, is_default*
Tabulated (mixed, conversion)   10    tab1, tab2, type, total_number
Tabulated scaled by lambda      11    filename*, max_force*
Tabulated (mixed, static)       12    tab1, tab2, mix value
Tabulated (cap radius)          13    filename, cap radius
Tabulated (scalled pairs)       14    filename, scale increment, max_force*
Lennard-Jones scaled by lambda  15    sigma*, epsilon*, max_force*
Lennard-Jones capped            16    sigma*, epsilon*, cap radius
==============================  ====  ======

Parameters with * are optional.
