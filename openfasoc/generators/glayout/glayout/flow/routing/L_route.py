from gdsfactory import cell
from gdsfactory.component import Component
from glayout.flow.routing.hashport import HPort
from glayout.flow.pdk.mappedpdk import MappedPDK
from typing import Optional, Union
from glayout.flow.primitives.via_gen import via_stack, via_array
from glayout.flow.pdk.util.comp_utils import evaluate_bbox, align_comp_to_port, to_decimal, move, to_float, prec_ref_center, get_primitive_rectangle, move
from glayout.flow.pdk.util.port_utils import rename_ports_by_orientation, rename_ports_by_list, print_ports, assert_port_manhattan, assert_ports_perpindicular
from decimal import Decimal
from gdsfactory.components.rectangle import rectangle as primitive_rectangle

def _get_port_coords(
	comp_bbox,
	comp_center,
	align_to: HPort,
	alignment: Optional[tuple[Optional[str],Optional[str]]] = None,
	port_width: Optional[float] = None,
	port_center: Optional[tuple[float,float]] = None
):
	"""Returns port coordinates for a component to align to another port"""
	# find center and bbox
	# import pdb; pdb.set_trace()
	cbbox = tuple(elem for elem in comp_bbox)
	if isinstance(comp_center, tuple):
		ccenter = comp_center 
	else:
		ccenter = tuple([comp_center.x, comp_center.y])
	# setup
	xdim = abs(cbbox[1][0] - cbbox[0][0])
	ydim = abs(cbbox[1][1] - cbbox[0][1])
	
	width = align_to.dwidth if not port_width else port_width
	alcenter = tuple(cc for cc in align_to.dcenter) if not port_center else port_center
	is_EW = bool(round(align_to.orientation + 90) % 180) # not EW == NS
	# error checks and decide orientation if None
	if alignment is None:
		if round(align_to.orientation) == 0:# facing east
			xalign = "r"
			yalign = "c"
		elif round(align_to.orientation) == 180:# facing west
			xalign = "l"
			yalign = "c"
		elif round(align_to.orientation) == 270:# facing south
			xalign = "c"
			yalign = "b"
		elif round(align_to.orientation) == 90:#facing north
			xalign = "c"
			yalign = "t"
		else:
			raise ValueError("port must be vertical or horizontal")
	else:
		xalign = (alignment[0] or "none").lower().strip()
		yalign = (alignment[1] or "none").lower().strip()
	# compute translation x amount for x alignment
	alignto_center = port_center if port_center else alcenter
	x_movcenter = alignto_center[0] - ccenter[0]
	x_mov_lr = abs(xdim/2 if is_EW else (width-xdim)/2)
	if "none" in xalign:
		xmov = 0
	elif "l" in xalign[0]:
		xmov = x_movcenter - x_mov_lr
	elif "c" in xalign[0]:
		xmov = x_movcenter
	elif "r" in xalign[0]:
		xmov = x_movcenter + x_mov_lr
	else:
		raise ValueError("please specify valid x alignment of l/r/c/None")
	# compute translation y amount for y alignment
	y_movcenter = alignto_center[1] - ccenter[1]
	y_move_updown = abs((width-ydim)/2 if is_EW else ydim/2)
	if "none" in yalign:
		ymov = 0
	elif "b" in yalign[0]:
		ymov = y_movcenter - y_move_updown
	elif "t" in yalign[0]:
		ymov = y_movcenter + y_move_updown
	elif "c" in yalign[0]:
		ymov = y_movcenter
	else:
		raise ValueError("please specify valid y alignment of t/b/c/None")
 
	retval = (xmov, ymov)
	return retval

# @cell(check_ports=False)
def L_route(
	pdk: MappedPDK,
	edge1: HPort,
	edge2: HPort,
	vwidth: Optional[float] = None,
	hwidth: Optional[float] = None,
	hglayer: Optional[str] = None,
	vglayer: Optional[str] = None,
	viaoffset: Optional[Union[tuple[bool,bool],bool]]=True,
	fullbottom: bool = True
) -> Component:
	"""creates a L shaped route between two Ports.
	
	edge1
	  |
	  ------|edge2
	
	REQUIRES: 
	- ports (a.k.a. edges) be vertical or horizontal
	- edges be perpindicular to each other
	
	DOES NOT REQUIRE:
	- correct 180 degree orientation of the port (e.g. a south facing port may result in north facing route)
	
	****NOTE: does no drc error checking (creates a dumb route)
	args:
	pdk = pdk to use
	edge1 = first port
	edge2 = second port
	vwidth = optional will default to vertical edge width if None
	hwidth = optional will default to horizontal edge width if None
	hglayer = glayer for vertical route. Defaults to the layer of the edge oriented N/S
	vglayer = glayer for horizontal route. Defaults to the layer of the edge oriented E/W
	viaoffset = push the via away from both edges so that inside corner aligns with via corner
	****via offset can also be specfied as a tuple(bool,bool): movex? if viaoffset[0] and movey? if viaoffset[1]
	fullbottom = fullbottom option for via
	"""
	# error checking, TODO: validate layers
	assert_port_manhattan([edge1,edge2])
	assert_ports_perpindicular(edge1,edge2)
	pdk.activate()
	Lroute = Component()
	# figure out which port is vertical
	vport = None
	hport = None
	edge1_is_EW = bool(round(edge1.orientation + 90) % 180)
	if edge1_is_EW:
		vport, hport = edge1, edge2
	else:
		hport, vport = edge1, edge2
	# arg setup
	vwidth = to_decimal(vwidth if vwidth else vport.dwidth)
	hwidth = to_decimal(hwidth if hwidth else hport.dwidth)
	hglayer = hglayer if hglayer else pdk.layer_to_glayer(vport.layer)
	vglayer = vglayer if vglayer else pdk.layer_to_glayer(hport.layer)
	if isinstance(viaoffset,bool):
		viaoffset = (True,True) if viaoffset else (False,False)
	# compute required dimensions
	hdim_center = to_decimal(vport.dcenter[0]) - to_decimal(hport.dcenter[0])
	vdim_center = to_decimal(hport.dcenter[1]) - to_decimal(vport.dcenter[1])
	hdim = abs(hdim_center) + hwidth/2
	vdim = abs(vdim_center) + vwidth/2
	# create and place vertical and horizontal connections
	# hconnect = Component()
	valign = ("l","c") if hdim_center > 0 else ("r","c")
	halign = ("c","b") if vdim_center > 0 else ("c","t")
	# vconnect = Component()
	temprect = get_primitive_rectangle(size=(hdim, vwidth))
	h_dest = _get_port_coords(temprect.bbox_np(), temprect.dcenter, vport, valign)
	temprect = get_primitive_rectangle(size=(hwidth, vdim))
	v_dest = _get_port_coords(temprect.bbox_np(), temprect.dcenter, hport, halign)
	hconnectref = get_primitive_rectangle(size=tuple(to_float((hdim,vwidth))), dest = h_dest,layer=pdk.get_glayer(hglayer))
	vconnectref = get_primitive_rectangle(size=tuple(to_float((hwidth,vdim))), dest = v_dest,layer=pdk.get_glayer(vglayer))
	# hconnect = get_primitive_rectangle(size=to_float((hdim,vwidth)),layer=pdk.get_glayer(hglayer))
	# vconnect = get_primitive_rectangle(size=to_float((hwidth,vdim)),layer=pdk.get_glayer(vglayer))
	#xalign
	#yalign
	
	# hconnectref = align_comp_to_port(hconnect, vport, valign, rtr_comp_ref=False)
	# hconnectref = move(hconnectref, (0, 0), (hx, hy))
	Lroute.add_ref(hconnectref)
	# vconnectref = align_comp_to_port(vconnect, hport, halign,  rtr_comp_ref=False)
	# vconnectref = move(vconnectref, (0, 0), (vx, vy))
	Lroute.add_ref(vconnectref)
	# create and place via (decide between via stack and via array)
	hv_via = Component()
	temp = hv_via << via_stack(pdk, hglayer, vglayer,fullbottom=fullbottom,fulltop=True)
	hv_via_dims = evaluate_bbox(hv_via,True)
	use_stack = hv_via_dims[0] > hwidth or hv_via_dims[1] > vwidth
	if not use_stack:
		temp = hv_via << via_array(pdk, hglayer, vglayer, size=tuple(to_float((hwidth,vwidth))), lay_bottom=True)
		hv_via.add_ports(temp.ports)
	hv_via_dims = evaluate_bbox(hv_via,True)
	h_to_v_via_ref = prec_ref_center(temp)
	h_to_v_via_ref.dmove(origin=(0,0),destination=(hport.dcenter[0], vport.dcenter[1]))
	Lroute.add(h_to_v_via_ref)
	if viaoffset[0] or viaoffset[1]:
		viadim_osx = hv_via_dims[0]/2
		viaxofs = abs(hwidth/2-viadim_osx)
		viaxofs = to_float(viaxofs if hdim_center > 0 else -1*viaxofs)
		viaxofs = viaxofs if viaoffset[0] else 0
		viadim_osy = hv_via_dims[1]/2
		viayofs = abs(vwidth/2-viadim_osy)
		viayofs = to_float(viayofs if vdim_center > 0 else -1*viayofs)
		viayofs = viayofs if viaoffset[1] else 0
		h_to_v_via_ref.dmovex(viaxofs).dmovey(viayofs)
	# add ports and return
	Lroute.add_ports(h_to_v_via_ref.ports)
	# Lroute.flatten()
	return rename_ports_by_orientation(Lroute)


