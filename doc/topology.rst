Topology file
=============

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