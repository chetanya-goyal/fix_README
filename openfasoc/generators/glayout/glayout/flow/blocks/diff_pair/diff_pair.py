from typing import Optional, Union

from gdsfactory import cell
from gdsfactory.component import Component, copy
from gdsfactory.components.rectangle import rectangle
from gdsfactory.routing.route_quad import route_quad
from gdsfactory.routing.route_sharp import route_sharp
from glayout.flow.pdk.mappedpdk import MappedPDK
from glayout.flow.routing.hashport import init_hashport
from glayout.flow.pdk.util.comp_utils import align_comp_to_port, evaluate_bbox, movex, movey
from glayout.flow.pdk.util.port_utils import (
    add_ports_perimeter,
    get_orientation,
    print_ports,
    rename_ports_by_list,
    rename_ports_by_orientation,
    set_port_orientation,
)
from glayout.flow.pdk.util.snap_to_grid import component_snap_to_grid
from glayout.flow.placement.common_centroid_ab_ba import common_centroid_ab_ba
from glayout.flow.primitives.fet import nmos, pmos
from glayout.flow.primitives.guardring import tapring
from glayout.flow.primitives.via_gen import via_stack
from glayout.flow.routing.c_route import c_route
from glayout.flow.routing.smart_route import smart_route
from glayout.flow.routing.straight_route import straight_route
from glayout.flow.spice import Netlist
from gdsfactory.snap import snap_to_grid


def diff_pair_netlist(fetL: Component, fetR: Component) -> Netlist:
	diff_pair_netlist = Netlist(circuit_name='DIFF_PAIR', nodes=['VP', 'VN', 'VDD1', 'VDD2', 'VTAIL', 'B'])
	diff_pair_netlist.connect_netlist(
		fetL.info['netlist'],
		[('D', 'VDD1'), ('G', 'VP'), ('S', 'VTAIL'), ('B', 'B')]
	)
	diff_pair_netlist.connect_netlist(
		fetR.info['netlist'],
		[('D', 'VDD2'), ('G', 'VN'), ('S', 'VTAIL'), ('B', 'B')]
	)
	return diff_pair_netlist

@cell(check_ports=False)
def diff_pair(
	pdk: MappedPDK,
	width: float = 3,
	fingers: int = 4,
	length: Optional[float] = None,
	n_or_p_fet: bool = True,
	plus_minus_seperation: float = 0,
	rmult: int = 1,
	dummy: Union[bool, tuple[bool, bool]] = True,
	substrate_tap: bool=True
) -> Component:
	"""create a diffpair with 2 transistors placed in two rows with common centroid place. Sources are shorted
	width = width of the transistors
	fingers = number of fingers in the transistors (must be 2 or more)
	length = length of the transistors, None or 0 means use min length
	short_source = if true connects source of both transistors
	n_or_p_fet = if true the diffpair is made of nfets else it is made of pfets
	substrate_tap: if true place a tapring around the diffpair (connects on met1)
	"""
	# TODO: error checking
	pdk.activate()
	diffpair = Component()
	# create transistors
	# import pdb;pdb.set_trace()
	well = None
	if isinstance(dummy, bool):
		dummy = (dummy, dummy)
	if n_or_p_fet:
		fetL = nmos(pdk, width=width, fingers=fingers,length=length,multipliers=1,with_tie=False,with_dummy=(dummy[0], False),with_dnwell=False,with_substrate_tap=False,rmult=rmult)
		fetR = nmos(pdk, width=width, fingers=fingers,length=length,multipliers=1,with_tie=False,with_dummy=(False,dummy[1]),with_dnwell=False,with_substrate_tap=False,rmult=rmult)
		min_spacing_x = pdk.get_grule("n+s/d")["min_separation"] - 2*(fetL.dxmax - fetL.ports["multiplier_0_plusdoped_E"].dcenter[0])
		well = "pwell"
	else:
		fetL = pmos(pdk, width=width, fingers=fingers,length=length,multipliers=1,with_tie=False,with_dummy=(dummy[0], False),dnwell=False,with_substrate_tap=False,rmult=rmult)
		fetR = pmos(pdk, width=width, fingers=fingers,length=length,multipliers=1,with_tie=False,with_dummy=(False,dummy[1]),dnwell=False,with_substrate_tap=False,rmult=rmult)
		min_spacing_x = pdk.get_grule("p+s/d")["min_separation"] - 2*(fetL.dxmax - fetL.ports["multiplier_0_plusdoped_E"].dcenter[0])
		well = "nwell"
	# place transistors
	viam2m3 = via_stack(pdk,"met2","met3",centered=True)
	metal_min_dim = max(pdk.get_grule("met2")["min_width"],pdk.get_grule("met3")["min_width"])
	metal_space = max(pdk.get_grule("met2")["min_separation"],pdk.get_grule("met3")["min_separation"],metal_min_dim)
	gate_route_os = evaluate_bbox(viam2m3)[0] - fetL.ports["multiplier_0_gate_W"].dwidth + metal_space
	min_spacing_y = metal_space + 2*gate_route_os
	min_spacing_y = min_spacing_y - 2*abs(fetL.ports["well_S"].dcenter[1] - fetL.ports["multiplier_0_gate_S"].dcenter[1])
	# TODO: fix spacing where you see +-0.5
	a_topl = (diffpair << fetL).dmovey(snap_to_grid(fetL.dymax+min_spacing_y/2+0.5, nm = 10)).dmovex(snap_to_grid(0-fetL.dxmax-min_spacing_x/2, nm = 10))
	b_topr = (diffpair << fetR).dmovey(snap_to_grid(fetR.dymax+min_spacing_y/2+0.5, nm = 10)).dmovex(snap_to_grid(fetL.dxmax+min_spacing_x/2, nm = 10))
	a_botr = (diffpair << fetR)
	a_botr.mirror_y().dmovey(snap_to_grid(0-0.5-fetL.dymax-min_spacing_y/2, nm = 10)).dmovex(snap_to_grid(fetL.dxmax+min_spacing_x/2, nm = 10))
	b_botl = (diffpair << fetL)
	b_botl.mirror_y().dmovey(snap_to_grid(0-0.5-fetR.dymax-min_spacing_y/2, nm = 10)).dmovex(snap_to_grid(0-fetL.dxmax-min_spacing_x/2, nm = 10))
	# if substrate tap place substrate tap
	if substrate_tap:
		diff_bbox = evaluate_bbox(diffpair, padding = 1)
		tapdim = tuple([snap_to_grid(diff_bbox[0], nm = 10), snap_to_grid(diff_bbox[1], nm = 10)])
		tapref = diffpair << tapring(pdk,tapdim,horizontal_glayer="met1")
		diffpair.add_ports(tapref.ports,prefix="tap_")
		try:
			port0 = init_hashport(a_topl.ports["multiplier_0_dummy_L_gsdcon_top_met_W"], layer_overload=pdk.get_glayer("met1"))
			port1 = init_hashport(diffpair.ports["tap_W_top_met_W"], layer_overload=pdk.get_glayer("met1"))
			diffpair<<straight_route(pdk,port0,port1,glayer2="met1")
		except KeyError:
			pass
		try:
			port0 = init_hashport(b_topr.ports["multiplier_0_dummy_R_gsdcon_top_met_W"], layer_overload=pdk.get_glayer("met1"))
			port1 = init_hashport(diffpair.ports["tap_E_top_met_E"], layer_overload=pdk.get_glayer("met1"))
			diffpair<<straight_route(pdk,port0,port1,glayer2="met1")
		except KeyError:
			pass
		try:
			port0 = init_hashport(b_botl.ports["multiplier_0_dummy_L_gsdcon_top_met_W"], layer_overload=pdk.get_glayer("met1"))
			port1 = init_hashport(diffpair.ports["tap_W_top_met_W"], layer_overload=pdk.get_glayer("met1"))
			diffpair<<straight_route(pdk,port0,port1,glayer2="met1")
		except KeyError:
			pass
		try:
			port0 = init_hashport(a_botr.ports["multiplier_0_dummy_R_gsdcon_top_met_W"], layer_overload=pdk.get_glayer("met1"))
			port1 = init_hashport(diffpair.ports["tap_E_top_met_E"], layer_overload=pdk.get_glayer("met1"))
			diffpair<<straight_route(pdk,port0,port1,glayer2="met1")
		except KeyError:
			pass
	# route sources (short sources)
	port0 = init_hashport(a_topl.ports["multiplier_0_source_E"], layer_overload=pdk.get_glayer("met1"))
	port1 = init_hashport(b_topr.ports["multiplier_0_source_W"], layer_overload=pdk.get_glayer("met2"))
	route_quad(diffpair, port0, port1, layer=pdk.get_glayer("met2"))
	port0 = init_hashport(b_botl.ports["multiplier_0_source_E"], layer_overload=pdk.get_glayer("met2"))
	port1 = init_hashport(a_botr.ports["multiplier_0_source_W"], layer_overload=pdk.get_glayer("met2"))
	route_quad(diffpair, port0, port1, layer=pdk.get_glayer("met2"))
	sextension = b_topr.ports["well_E"].dcenter[0] - b_topr.ports["multiplier_0_source_E"].dcenter[0]
	port0 = init_hashport(b_topr.ports["multiplier_0_source_E"], layer_overload=pdk.get_glayer("met2"))
	port1 = init_hashport(a_botr.ports["multiplier_0_source_E"], layer_overload=pdk.get_glayer("met2"))
	source_routeE = diffpair << c_route(pdk, port0, port1,extension=sextension)
	port0 = init_hashport(a_topl.ports["multiplier_0_source_W"], layer_overload=pdk.get_glayer("met2"))
	port1 = init_hashport(b_botl.ports["multiplier_0_source_W"], layer_overload=pdk.get_glayer("met2"))
	source_routeW = diffpair << c_route(pdk, port0, port1,extension=sextension)
	# route drains
	# place via at the drain
	drain_br_via = diffpair << viam2m3
	drain_bl_via = diffpair << viam2m3
	center_snapr = tuple([snap_to_grid(a_botr.ports["multiplier_0_drain_N"].dcenter[0], nm = 10), snap_to_grid(a_botr.ports["multiplier_0_drain_N"].dcenter[1], nm = 10)])
	center_snapl = tuple([snap_to_grid(b_botl.ports["multiplier_0_drain_N"].dcenter[0], nm = 10), snap_to_grid(b_botl.ports["multiplier_0_drain_N"].dcenter[1], nm = 10)])
	via_snapped_ydim = snap_to_grid(viam2m3.dymin, nm = 10)
	drain_br_via.dmove(center_snapr).dmovey(via_snapped_ydim)
	drain_bl_via.dmove(center_snapl).dmovey(via_snapped_ydim)
	drain_br_viatm = diffpair << viam2m3
	drain_bl_viatm = diffpair << viam2m3
	drain_br_viatm.dmove(center_snapr).dmovey(via_snapped_ydim)
	drain_bl_viatm.dmove(center_snapl).dmovey(snap_to_grid(-1.5 * evaluate_bbox(viam2m3)[1] - metal_space, nm = 10))
	# create route to drain via
	width_drain_route = b_topr.ports["multiplier_0_drain_E"].dwidth
	dextension = source_routeE.dxmax - b_topr.ports["multiplier_0_drain_E"].dcenter[0] + metal_space
	bottom_extension = viam2m3.dymax + width_drain_route/2 + 2*metal_space
	drain_br_viatm.dmovey(snap_to_grid(0-bottom_extension - metal_space - width_drain_route/2 - viam2m3.dymax, nm = 10))
	port0 = init_hashport(drain_br_viatm.ports["top_met_N"], layer_overload=pdk.get_glayer("met3"))
	port1 = init_hashport(drain_br_via.ports["top_met_S"], layer_overload=pdk.get_glayer("met3"))
	route_quad(diffpair, port0, port1, layer=pdk.get_glayer("met3"))
	port0 = init_hashport(drain_bl_viatm.ports["top_met_N"], layer_overload=pdk.get_glayer("met3"))
	port1 = init_hashport(drain_bl_via.ports["top_met_S"], layer_overload=pdk.get_glayer("met3"))
	route_quad(diffpair, port0, port1, layer=pdk.get_glayer("met3"))
	p0 = init_hashport(drain_bl_via.ports["bottom_met_W"], layer_overload=pdk.get_glayer("met2")).move_copy(drain_bl_via.ports['bottom_met_W'].dcenter[0], drain_bl_via.ports['bottom_met_W'].dcenter[1] - bottom_extension)
	p1 = init_hashport(drain_br_via.ports["bottom_met_E"], layer_overload=pdk.get_glayer('met2')).move_copy(drain_br_via.ports['bottom_met_E'].dcenter[0], drain_br_via.ports['bottom_met_E'].dcenter[1] - bottom_extension - metal_space - width_drain_route)
	floating_port_drain_bottom_L = set_port_orientation(p0, get_orientation("E"))
	floating_port_drain_bottom_R = set_port_orientation(p1, get_orientation("W"))
	port1 = init_hashport(b_topr.ports["multiplier_0_drain_E"], layer_overload=pdk.get_glayer('met2'))
	drain_routeTR_BL = diffpair << c_route(pdk, floating_port_drain_bottom_L, port1,extension=dextension, width1=width_drain_route,width2=width_drain_route)
	port1 = init_hashport(a_topl.ports["multiplier_0_drain_W"], layer_overload=pdk.get_glayer('met2'))
	drain_routeTL_BR = diffpair << c_route(pdk, floating_port_drain_bottom_R, port1,extension=dextension, width1=width_drain_route,width2=width_drain_route)
	# cross gate route top with c_route. bar_minus ABOVE bar_plus
	get_left_extension = lambda bar, a_topl=a_topl, diffpair=diffpair, pdk=pdk : (abs(diffpair.dxmin-min(a_topl.ports["multiplier_0_gate_W"].dcenter[0],bar.ports["e1"].dcenter[0])) + pdk.get_grule("met2")["min_separation"])
	get_right_extension = lambda bar, b_topr=b_topr, diffpair=diffpair, pdk=pdk : (abs(diffpair.dxmax-max(b_topr.ports["multiplier_0_gate_E"].dcenter[0],bar.ports["e3"].dcenter[0])) + pdk.get_grule("met2")["min_separation"])
	# lay bar plus and PLUSgate_routeW
	bar_comp = rectangle(centered=True,size=(abs(snap_to_grid(b_topr.dxmax-a_topl.dxmin, nm = 10)), snap_to_grid(b_topr.ports["multiplier_0_gate_E"].dwidth, nm = 10)),layer=pdk.get_glayer("met2"))
	bar_plus = (diffpair << bar_comp).dmovey(snap_to_grid(diffpair.dymax + bar_comp.dymax + pdk.get_grule("met2")["min_separation"], nm = 10))
	PLUSgate_routeW = diffpair << c_route(pdk, init_hashport(a_topl.ports["multiplier_0_gate_W"], layer_overload=pdk.get_glayer('met2')), init_hashport(bar_plus.ports["e1"], layer_overload=pdk.get_glayer('met2')), extension=get_left_extension(bar_plus))
	# lay bar minus and MINUSgate_routeE
	plus_minus_seperation = max(pdk.get_grule("met2")["min_separation"], plus_minus_seperation)
	bar_minus = (diffpair << bar_comp).dmovey(snap_to_grid(diffpair.dymax +bar_comp.dymax + plus_minus_seperation, nm = 10))
	MINUSgate_routeE = diffpair << c_route(pdk, init_hashport(b_topr.ports["multiplier_0_gate_E"], layer_overload=pdk.get_glayer('met2')), init_hashport(bar_minus.ports["e3"], layer_overload=pdk.get_glayer('met2')), extension=get_right_extension(bar_minus))
	# lay MINUSgate_routeW and PLUSgate_routeE
	MINUSgate_routeW = diffpair << c_route(pdk, set_port_orientation(init_hashport(b_botl.ports["multiplier_0_gate_E"], layer_overload=pdk.get_glayer('met2')),"W"), init_hashport(bar_minus.ports["e1"], layer_overload=pdk.get_glayer('met2')), extension=get_left_extension(bar_minus))
	PLUSgate_routeE = diffpair << c_route(pdk, set_port_orientation(init_hashport(a_botr.ports["multiplier_0_gate_W"], layer_overload=pdk.get_glayer('met2')),"E"), init_hashport(bar_plus.ports["e3"], layer_overload=pdk.get_glayer('met2')), extension=get_right_extension(bar_plus))
	# correct pwell place, add ports, flatten, and return
	diffpair.add_ports(a_topl.ports,prefix="tl_")
	diffpair.add_ports(b_topr.ports,prefix="tr_")
	diffpair.add_ports(b_botl.ports,prefix="bl_")
	diffpair.add_ports(a_botr.ports,prefix="br_")
	diffpair.add_ports(source_routeE.ports,prefix="source_routeE_")
	diffpair.add_ports(source_routeW.ports,prefix="source_routeW_")
	diffpair.add_ports(drain_routeTR_BL.ports,prefix="drain_routeTR_BL_")
	diffpair.add_ports(drain_routeTL_BR.ports,prefix="drain_routeTL_BR_")
	diffpair.add_ports(MINUSgate_routeW.ports,prefix="MINUSgateroute_W_")
	diffpair.add_ports(MINUSgate_routeE.ports,prefix="MINUSgateroute_E_")
	diffpair.add_ports(PLUSgate_routeW.ports,prefix="PLUSgateroute_W_")
	diffpair.add_ports(PLUSgate_routeE.ports,prefix="PLUSgateroute_E_")
	# diffpair.add_padding(layers=(pdk.get_glayer(well),), default=0)
	component = component_snap_to_grid(rename_ports_by_orientation(diffpair))

	# component.info['netlist'] = diff_pair_netlist(fetL, fetR)
	return component



@cell
def diff_pair_generic(
	pdk: MappedPDK,
	width: float = 3,
	fingers: int = 4,
	length: Optional[float] = None,
	n_or_p_fet: bool = True,
	plus_minus_seperation: float = 0,
	rmult: int = 1,
	dummy: Union[bool, tuple[bool, bool]] = True,
	substrate_tap: bool=True
) -> Component:
	diffpair = common_centroid_ab_ba(pdk,width,fingers,length,n_or_p_fet,rmult,dummy,substrate_tap)
	diffpair << smart_route(pdk,diffpair.ports["A_source_E"],diffpair.ports["B_source_E"],diffpair, diffpair)
	return diffpair
