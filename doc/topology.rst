Functionals
-----------

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


Topology file
-------------

[ bondtypes ]
+++++++++++++

========================  =====  =======
Name of interaction       func   params
========================  =====  =======
Harmonic eq1_             1      K [1]_, r0
FENE eqFENE_              7      K, r0
Tabulated                 8      table index
========================  =====  =======

.. [1] Force constant internally divided by 2.0

[ angletypes ]
++++++++++++++

========================  =====  =======
Name of interaction       func   params
========================  =====  =======
Harmonic eq2_             1      K [2]_, theta0 (deg)
Tabulated                 8      table index
Cosine   eq3_             11     K [2]_, theta0 (deg)
========================  =====  =======

.. [2] Force constant internally divided by 2.0

[ dihedraltypes ]
+++++++++++++++++

========================  =====  =======
Name of interaction       func   params
========================  =====  =======
HarmonicNCos  eq4_        1      K, phi0 (deg), multiplicity
Ryckaert Bellemans  eq5_  3      K0, K1, K2, K3, K4, K5
Tabulated                 8      table index
========================  =====  =======


[ nonbond_params ]
++++++++++++++++++

==============================  ====  ======
Name of interaction             func  params
==============================  ====  ======
Lennard-Jones                   1     sigma*, epsilon*
Tabulated                       8     filename*
Tabulated (conversion)          9     filename*, type, total number, p_min, p_max, is_default*
Tabulated (mixed, conversion)   10    tab1, tab2, type, total_number
Tabulated scaled by lambda      11    filename*, max_force*
Tabulated (mixed)               12    tab1, tab2, mix value
Tabulated (cap radius)          13    filename, cap radius
Tabulated (scalled pairs)       14    filename, scale increment, max_force*
Lennard-Jones scaled by lambda  15    sigma*, epsilon*, max_force*
Lennard-Jones capped            16    sigma*, epsilon*, cap radius
==============================  ====  ======

Parameters with * are optional.
