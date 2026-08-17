
Ring Modulator Layout Program 2
===============================================

Marc Schneider, 2026

marc.schneider@kit.edu

This program is for designing silicon photonic ring modulators (RM) with a
pn-junction along the circular ring to modulate light traversing through a
loosely coupled, straight bus waveguide. The result can then be saved as GDSII
file for use in your photonic chip design.
In resonance of the ring, light can be coupled to a straight drop waveguide,
whose distance to the ring and therefore coupling parameter can be adjusted
independently from the distance of the bus waveguide to the ring. The outer
waveguide connections are designed as strip waveguides, which taper to rib
waveguides with slabs in the whole RM region. The slab thickness of the bus
and drop waveguides, as well as the outer ring slab is defined by the GDS
layer for the dedicated etch depth. The slab thickness of the inner ring slab
can be chosen as the same or a different GDS layer with a different etch depth.
To change the resonance of the ring through the thermo-optic effect, a heater
on top of the ring is implemented.
Some geometry aspects are influenced by the design rules of a certain European
fab, but as the design is widely adjustable, it should be usable for many
different fabs.
The program is written in Python (my first Python program, so you might notice
some awkward constructions) and depends heavily on the gdsfactory library
(https://github.com/gdsfactory/gdsfactory).


The project is published under the MIT license (see LICENSE), with the
following exceptions:

- The included logo pictures are not covered by the MIT license.

This software was developed at the Karlsruhe Institute of Technology (KIT),
Germany. This software is an experimental system. KIT or the author assume no
responsibility whatsoever for its use by other parties, and makes no
guarantees, expressed or implied, about its quality, reliability, or any
other characteristics.




Help
----
If you find a bug in the software or need assistance, please send a message to 
marc.schneider@kit.edu.
