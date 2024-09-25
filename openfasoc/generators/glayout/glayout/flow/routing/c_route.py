from gdsfactory import cell
from gdsfactory.component import Component
from glayout.flow.routing.hashport import HPort, init_hashport
from glayout.flow.pdk.mappedpdk import MappedPDK
from typing import Optional, Union
from math import isclose
from glayout.flow.primitives.via_gen import via_stack
from glayout.flow.routing.straight_route import straight_route
from gdsfactory.components.rectangle import rectangle
from glayout.flow.pdk.util.comp_utils import evaluate_bbox, get_primitive_rectangle
from glayout.flow.pdk.util.port_utils import add_ports_perimeter, rename_ports_by_orientation, rename_ports_by_list, print_ports, set_port_width, set_port_orientation, get_orientation
from pydantic import validate_arguments
from gdsfactory.snap import snap_to_grid
from gdsfactory import get_layer

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def __fill_empty_viastack__macro(pdk: MappedPDK, glayer: str, size: Optional[tuple[float,float]]=None) -> Component:
	"""returns a rectangle with ports that pretend to be viastack ports
	by default creates a rectangle with size double the min width of the glayer"""
	if size is None:
		size = pdk.snap_to_2xgrid([2*pdk.get_grule(glayer)["min_width"],2*pdk.get_grule(glayer)["min_width"]])
	comp = rectangle(size=size,layer=pdk.get_glayer(glayer),centered=True)
	comp = rename_ports_by_orientation(rename_ports_by_list(comp,replace_list=[("e","top_met_")]))
	return comp

@cell(check_ports=False)
def c_route(
	pdk: MappedPDK, 
	edge1: HPort, 
	edge2: HPort, 
	extension: Optional[float]=0.5, 
	width1: Optional[float] = None, 
	width2: Optional[float] = None,
	cwidth: Optional[float] = None,
	e1glayer: Optional[str] = None, 
	e2glayer: Optional[str] = None, 
	cglayer: Optional[str] = None, 
	viaoffset: Optional[Union[bool,tuple[Optional[bool],Optional[bool]]]]=(True,True),
	fullbottom: Optional[bool] = False,
	debug=False
) -> Component:
	"""creates a C shaped route between two Ports.
	
	edge1--|
	       |
	edge2--|
	
	REQUIRES: ports be parralel vertical or horizontal edges
	****NOTE: does no drc error checking (creates a dumb route)
	args:
	pdk = pdk to use
	edge1 = first port
	edge2 = second port
	width1 = optional will default to edge1 width if None
	width2 = optional will default to edge2 width if None
	e1glayer = glayer for the parts connecting to the edge1. Default to layer of edge1
	e2glayer = glayer for the parts connecting to the edge2. Default to layer of edge2
	cglayer = glayer for the connection part (part that goes through a via) defaults to e1glayer met+1
	viaoffset = offsets the via so that it is flush with the cglayer (may be needed for drc) i.e. -| vs _|
	- True offsets via towards the other via
	- False offsets via away from the other via
	- None means center (no offset)
	****NOTE: viaoffset pushes both vias towards each other slightly
	"""
	if debug:
		pass
	# error checking and figure out args
	if round(edge1.orientation) % 90 or round(edge2.orientation) % 90:
		raise ValueError("Ports must be vertical or horizontal")
	if not isclose(edge1.orientation,edge2.orientation):
		raise ValueError("Ports must be parralel and have same orientation")
	width1 = width1 if width1 else edge1.dwidth
	width2 = width2 if width2 else edge1.dwidth
	e1glayer = e1glayer if e1glayer else pdk.layer_to_glayer(edge1.layer)
	e2glayer = e2glayer if e2glayer else pdk.layer_to_glayer(edge2.layer)
	eglayer_plusone = "met" + str(int(e1glayer[-1])+1)
	cglayer = cglayer if cglayer else eglayer_plusone
	if not "met" in e1glayer or not "met" in e2glayer or not "met" in cglayer:
		raise ValueError("given layers must be metals")
	viaoffset = (None, None) if viaoffset is None else viaoffset
	if isinstance(viaoffset,bool):
		viaoffset = (True,True) if viaoffset else (False,False)
	pdk.has_required_glayers([e1glayer,e2glayer,cglayer])
	pdk.activate()
	extension = snap_to_grid(extension)
	# create route
	croute = Component()
	viastack1 = via_stack(pdk,e1glayer,cglayer,fullbottom=fullbottom,assume_bottom_via=True)
	viastack2 = via_stack(pdk,e2glayer,cglayer,fullbottom=fullbottom,assume_bottom_via=True)
	if e1glayer==cglayer and e2glayer==cglayer:
		viastack1 = __fill_empty_viastack__macro(pdk,e1glayer)
		viastack2 = __fill_empty_viastack__macro(pdk,e2glayer)
	elif e1glayer == cglayer:
		viastack1 = __fill_empty_viastack__macro(pdk,e1glayer,size=evaluate_bbox(viastack2))
	elif e2glayer == cglayer:
		viastack2 = __fill_empty_viastack__macro(pdk,e2glayer,size=evaluate_bbox(viastack1))
	# find extension
	e1_length = snap_to_grid(extension + evaluate_bbox(viastack1)[0])
	e2_length = snap_to_grid(extension + evaluate_bbox(viastack2)[0])
	xdiff = snap_to_grid(abs(edge1.dcenter[0] - edge2.dcenter[0]))
	ydiff = snap_to_grid(abs(edge1.dcenter[1] - edge2.dcenter[1]))
	if not isclose(edge1.dcenter[0],edge2.dcenter[0]):
		if round(edge1.orientation) == 0:# facing east
			if edge1.dcenter[0] > edge2.dcenter[0]:
				e2_length += xdiff
			else:
				e1_length += xdiff
		elif round(edge1.orientation) == 180:# facing west
			if edge1.dcenter[0] < edge2.dcenter[0]:
				e2_length += xdiff
			else:
				e1_length += xdiff
	if not isclose(edge1.dcenter[1],edge2.dcenter[1]):
		if round(edge1.orientation) == 270:# facing south
			if edge1.dcenter[1] < edge2.dcenter[1]:
				e2_length += ydiff
			else:
				e1_length += ydiff
		elif round(edge1.orientation) == 90:#facing north
			if edge1.dcenter[1] > edge2.dcenter[1]:
				e2_length += ydiff
			else:
				e1_length += ydiff
	# move into position
	e1_extension_comp = Component()
	e2_extension_comp = Component()
	box_dims = [(e1_length, width1),(e2_length, width2)]
	if round(edge1.orientation) == 90 or round(edge1.orientation) == 270:
		box_dims = [(width1, e1_length),(width2, e2_length)]
	rect_c1 = get_primitive_rectangle(size=box_dims[0], layer=pdk.get_glayer(e1glayer))
	rect_c2 = get_primitive_rectangle(size=box_dims[1], layer=pdk.get_glayer(e2glayer))
	rect_c1 = rename_ports_by_orientation(rename_ports_by_list(rect_c1,[("e","e_")]))
	rect_c2 = rename_ports_by_orientation(rename_ports_by_list(rect_c2,[("e","e_")]))
	# TODO: make sure ports match bbox
	e1_extension = e1_extension_comp << rect_c1
	e2_extension = e2_extension_comp << rect_c2
	e1_extension.dmove(origin = (0,0), destination=edge1.dcenter)
	e2_extension.dmove(origin = (0,0), destination=edge2.dcenter)
	if round(edge1.orientation) == 0:# facing east
		e1_extension.dmovey(origin = 0, destination = 0-evaluate_bbox(e1_extension)[1]/2)
		e2_extension.dmovey(origin = 0, destination = 0-evaluate_bbox(e2_extension)[1]/2)
	elif round(edge1.orientation) == 180:# facing west
		e1_extension.dmovex(origin = 0, destination = 0-evaluate_bbox(e1_extension)[0])
		e2_extension.dmovex(origin = 0, destination = 0-evaluate_bbox(e2_extension)[0])
		e1_extension.dmovey(origin = 0, destination = 0-evaluate_bbox(e1_extension)[1]/2)
		e2_extension.dmovey(origin = 0, destination = 0-evaluate_bbox(e2_extension)[1]/2)
	elif round(edge1.orientation) == 270:# facing south
		e1_extension.dmovex(origin = 0, destination = 0-evaluate_bbox(e1_extension)[0]/2)
		e2_extension.dmovex(origin = 0, destination = 0-evaluate_bbox(e2_extension)[0]/2)
		e1_extension.dmovey(origin = 0, destination = 0-evaluate_bbox(e1_extension)[1])
		e2_extension.dmovey(origin = 0, destination = 0-evaluate_bbox(e2_extension)[1])
	else:#facing north
		e1_extension.dmovex(origin = 0, destination = 0-evaluate_bbox(e1_extension)[0]/2)
		e2_extension.dmovex(origin = 0, destination = 0-evaluate_bbox(e2_extension)[0]/2)
	# place viastacks
	e1_extension_comp.add_ports(e1_extension.ports)
	e2_extension_comp.add_ports(e2_extension.ports)
	me1 = e1_extension_comp << viastack1
	me2 = e2_extension_comp << viastack2
	route_ports = [None,None]
	via_flush = snap_to_grid(abs((width1 - evaluate_bbox(viastack1)[0])/2) if viaoffset else 0)
	via_flush1 = via_flush if viaoffset[0] else 0-via_flush
	via_flush1 = 0 if viaoffset[0] is None else via_flush1
	via_flush2 = via_flush if viaoffset[1] else 0-via_flush
	via_flush2 = 0 if viaoffset[1] is None else via_flush2
	if round(edge1.orientation) == 0:# facing east
		me1.dmove(origin = (0,0), destination=e1_extension.ports["e_E"].dcenter)
		me2.dmove(origin = (0,0), destination=e2_extension.ports["e_E"].dcenter)
		via_flush *= 1 if me1.dymax > me2.dymax else -1
		me1.dmovex(origin = 0, destination = 0-viastack1.dxmax).dmovey(origin = 0, destination = 0-via_flush1)
		me2.dmovex(origin = 0, destination = 0-viastack2.dxmax).dmovey(origin = 0, destination = via_flush2)
		me1, me2 = (me1, me2) if (me1.dcenter.y > me2.dcenter.y) else (me2, me1)
		route_ports = [init_hashport(me1.ports["top_met_N"], layer_overload = pdk.get_glayer('met3')),init_hashport(me2.ports["top_met_S"], layer_overload = pdk.get_glayer('met3'))]
		fix_connection_direction = "E"
		fix_ports = [init_hashport(me1.ports["top_met_E"], layer_overload = pdk.get_glayer('met3')),init_hashport(me2.ports["top_met_E"], layer_overload = pdk.get_glayer('met3'))]
	elif round(edge1.orientation) == 180:# facing west
		me1.dmove(origin = (0,0), destination=e1_extension.ports["e_W"].dcenter)
		me2.dmove(origin = (0,0), destination=e2_extension.ports["e_W"].dcenter)
		via_flush *= 1 if me1.dymax > me2.dymax else -1
		me1.dmovex(origin = 0, destination = viastack1.dxmax).dmovey(origin = 0, destination = 0-via_flush1)
		me2.dmovex(origin = 0, destination = viastack2.dxmax).dmovey(origin = 0, destination = via_flush2)
		me1, me2 = (me1, me2) if (me1.dcenter.y > me2.dcenter.y) else (me2, me1)
		route_ports = [init_hashport(me1.ports["top_met_N"], layer_overload = pdk.get_glayer('met3')),init_hashport(me2.ports["top_met_S"], layer_overload = pdk.get_glayer('met3'))]
		fix_connection_direction = "E"
		fix_ports = [init_hashport(me1.ports["top_met_E"], layer_overload = pdk.get_glayer('met3')),init_hashport(me2.ports["top_met_E"], layer_overload = pdk.get_glayer('met3'))]
	elif round(edge1.orientation) == 270:# facing south
		me1.dmove(origin = (0,0), destination=e1_extension.ports["e_S"].dcenter)
		me2.dmove(origin = (0,0), destination=e2_extension.ports["e_S"].dcenter)
		via_flush *= 1 if me1.dxmax > me2.dxmax else -1
		me1.dmovey(origin = 0, destination = viastack1.dxmax).dmovex(origin = 0, destination = 0-via_flush1)
		me2.dmovey(origin = 0, destination = viastack2.dxmax).dmovex(origin = 0, destination = via_flush2)
		me1, me2 = (me1, me2) if (me1.dcenter.x > me2.dcenter.x) else (me2, me1)
		route_ports = [init_hashport(me1.ports["top_met_E"], layer_overload = pdk.get_glayer('met3')),init_hashport(me2.ports["top_met_W"], layer_overload = pdk.get_glayer('met3'))]
		fix_connection_direction = "N"
		fix_ports = [init_hashport(me1.ports["top_met_N"], layer_overload = pdk.get_glayer('met3')),init_hashport(me2.ports["top_met_N"], layer_overload = pdk.get_glayer('met3'))]
	else:#facing north
		me1.dmove(origin = (0,0), destination=e1_extension.ports["e_N"].dcenter)
		me2.dmove(origin = (0,0), destination=e2_extension.ports["e_N"].dcenter)
		via_flush *= 1 if me1.dxmax > me2.dxmax else -1
		me1.dmovey(origin = 0, destination = 0-viastack1.dxmax).dmovex(origin = 0, destination = 0-via_flush1)
		me2.dmovey(origin = 0, destination = 0-viastack2.dxmax).dmovex(origin = 0, destination = via_flush2)
		me1, me2 = (me1, me2) if (me1.dcenter.x > me2.dcenter.x) else (me2, me1)
		route_ports = [init_hashport(me1.ports["top_met_E"], layer_overload = pdk.get_glayer('met3')),init_hashport(me2.ports["top_met_W"], layer_overload = pdk.get_glayer('met3'))]
		fix_connection_direction = "N"
		fix_ports = [init_hashport(me1.ports["top_met_N"], layer_overload = pdk.get_glayer('met3')),init_hashport(me2.ports["top_met_N"], layer_overload = pdk.get_glayer('met3'))]
	# connect extensions, add ports, return
	croute << e1_extension_comp
	croute << e2_extension_comp
	if cwidth:
		route_ports = [set_port_width(port_,cwidth) for port_ in route_ports]
	maxwidth = max(route_ports[0].width, route_ports[1].width)
	# route_ports[0].dwidth = route_ports[1].dwidth = max(route_ports[0].dwidth, route_ports[1].dwidth)
	# import pdb; pdb.set_trace()
	route_port0 = init_hashport(route_ports[0], width_overload=maxwidth, layer_overload=pdk.get_glayer(cglayer))
	route_port1 = init_hashport(route_ports[1], width_overload=maxwidth, layer_overload=pdk.get_glayer(cglayer))
	cconnection = croute << straight_route(pdk, route_port0,route_port1,glayer1=cglayer,glayer2=cglayer)
	for _port in fix_ports:
		port2 = init_hashport(cconnection.ports["route_"+fix_connection_direction], layer_overload=pdk.get_glayer(cglayer))
		port1 = init_hashport(_port, layer_overload=pdk.get_glayer(cglayer))
		# _port.layer = pdk.get_glayer(cglayer)
		# print("\n\n")
		# print(_port.name, port2.name)
		# print("\n\n")
		croute << straight_route(pdk, port1, port2, glayer1=cglayer,glayer2=cglayer)
	new_port_table = []
	orta1 = get_orientation(route_ports[0].orientation, int_only=True)
	orta2 = get_orientation(route_ports[1].orientation, int_only=True)
	new_port_table = [init_hashport(route_ports[0], orientation_overload=float(orta1), layer_overload=get_layer(pdk.get_glayer(cglayer)), width_overload=maxwidth),init_hashport(route_ports[1], orientation_overload=orta2, layer_overload=get_layer(pdk.get_glayer(cglayer)), width_overload=maxwidth)]
	# for i,port_to_add in enumerate(route_ports):
	# 	orta = get_orientation(port_to_add.orientation)
	# 	route_ports[i] = set_port_orientation(port_to_add, orta)
	croute.add_ports(new_port_table,prefix="con_")
	# croute.flatten()
	return rename_ports_by_orientation(rename_ports_by_list(croute, [("con_","con_")]))

