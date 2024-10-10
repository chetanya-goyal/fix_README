from glayout.flow.pdk.mappedpdk import MappedPDK
from gdsfactory import cell
from gdsfactory.component import Component
from gdsfactory.components.rectangle import rectangle
# from gdsfactory.components.rectangular_ring import rectangular_ring
from glayout.flow.pdk.util.comp_utils import create_rectangular_ring
from glayout.flow.primitives.via_gen import via_array, via_stack
from typing import Optional
from glayout.flow.pdk.util.comp_utils import to_decimal, to_float, evaluate_bbox
from glayout.flow.pdk.util.port_utils import print_ports
from glayout.flow.pdk.util.snap_to_grid import component_snap_to_grid
from glayout.flow.routing.L_route import L_route
from glayout.flow.primitives.mimcap import mimcap, mimcap_array
from glayout.flow.primitives.fet import nmos, pmos, multiplier
from glayout.flow.pdk.sky130_mapped import sky130_mapped_pdk as sky130
from glayout.flow.primitives.guardring import tapring
from glayout.flow.blocks.diff_pair import diff_pair
from glayout.flow.placement.two_transistor_place import two_transistor_place
from glayout.flow.placement.two_transistor_interdigitized import two_transistor_interdigitized, two_pfet_interdigitized, two_nfet_interdigitized
from glayout.flow.placement.four_transistor_interdigitized import generic_4T_interdigitzed 
from glayout.flow.blocks.current_mirror import current_mirror

# celll = nmos(sky130, multipliers=2)
# celll = diff_pair(sky130)
# celll = multiplier(sky130, "n+s/d")
# celll = via_stack(sky130, 'active_tap', 'met2', fullbottom=True)
# celll = tapring(sky130, (5, 5), sides=(True, True, True, True))
# celll = mimcap(sky130, (5, 5))
# celll = mimcap_array(sky130, 5, 5, (5, 5))
# celll.write_gds()
# celll = via_array(sky130, 'met1', 'met2', (5, 5))
# celll = via_stack(sky130, 'met1', 'met2')
# celll = two_transistor_interdigitized(sky130, 'pfet', 5)
# celll = two_transistor_place(sky130, 'aba bab', (nmos, {'pdk': sky13 0}), (pmos, {'pdk': sky130}))
# celll = generic_4T_interdigitzed(sky130, 'nfet', 'pfet', 5)
# celll = current_mirror(sky130, 5, 'nfet')

# copied = cell.copy()
# copied.name = 'idk'
celll.write_gds('holyshit.gds')

# sky130.drc_magic(copied, 'idk')