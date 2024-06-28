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

cell = nmos(sky130)
# cell = multiplier(sky130, "n+s/d")
# cell = via_stack(sky130, 'active_tap', 'met2', fullbottom=True)
# cell = tapring(sky130, (5, 5), sides=(True, True, True, True))
# cell = mimcap(sky130, (5, 5))
# cell = mimcap_array(sky130, 5, 5, (5, 5))
# cell = via_stack(sky130, 'met1', 'met2')
# cell = via_array(sky130, 'met1', 'met2', (5, 5))

copied = cell.copy()
copied.name = 'idk'
cell.write_gds('output.gds')

sky130.drc_magic(copied, 'idk')