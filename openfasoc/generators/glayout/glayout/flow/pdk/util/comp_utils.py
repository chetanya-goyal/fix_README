from pydantic import validate_arguments
from gdsfactory import Component, ComponentReference
from gdsfactory.snap import snap_to_grid
from gdsfactory.components import bbox
# from gdsfactory.typings import Component, ComponentReference
from gdsfactory.components.rectangle import rectangle
from glayout.flow.routing.hashport import HPort
from gdsfactory.components import rectangle

from typing import Callable, Union, Optional, Iterable
from decimal import Decimal
# from gdsfactory.functions import transformed
# from gdsfactory.functions import move as __gf_move
from glayout.flow.pdk.mappedpdk import MappedPDK
from gdstk import rectangle as primitive_rectangle
from .port_utils import add_ports_perimeter, rename_ports_by_list, parse_direction
from gdsfactory.typings import LayerSpec
import gdsfactory.boolean as boolean
from gdsfactory.snap import snap_to_grid2x
import kfactory as kf
import gdsfactory as gf
def transformed(ref: ComponentReference):
    """Returns flattened cell with reference transformations applied.

    Args:
        ref: the reference to flatten into a new cell.

    """
    # from gdsfactory.component import copy_reference
    # import pdb; pdb.set_trace()
    c = Component()
    c = ref.cell.from_kcell()
    
    c.flatten()
    # c.copy_child_info(ref.parent_cell)
    c.add_ports(ref.ports)
    c.info["transformed_cell"] = ref.parent_cell.name
    return c

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def evaluate_bbox(custom_comp: Union[Component, ComponentReference], return_decimal: Optional[bool]=False, padding: float=0) -> tuple[Union[float,Decimal],Union[float,Decimal]]:
	"""returns the length and height of a component like object"""
	import kfactory as kf
	if isinstance(custom_comp, ComponentReference) or isinstance(custom_comp, kf.Instance):
		newcomp = transformed(custom_comp)
		compbbox = newcomp.bbox_np()
	else:
		compbbox = custom_comp.bbox_np()
	width = abs(Decimal(str(compbbox[1][0])) - Decimal(str(compbbox[0][0]))) + 2*Decimal(str(padding))
	height = abs(Decimal(str(compbbox[1][1])) - Decimal(str(compbbox[0][1]))) + 2*Decimal(str(padding))
	if return_decimal:
		return (width,height)
	return (float(width),float(height))
	

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def center_to_edge_distance(custom_comp: Union[Component, ComponentReference], direction: Union[str,int]) -> float:
	"""specifies the distance between the center of custom_comp and a specified edge given by the direction argument
	the component is considerd a rectangle (using the bounding box), such there are 4 edges
	Args:
		custom_comp (Component | ComponentReference): Component we want the center of
		direction (str | int): the edge we are interested in
	Returns:
		float: absolute distance between custom_comp center and N,S,E, or W edge
	"""
	compbbox = custom_comp.bbox_np()
	#center = prec_center(custom_comp)
	center = custom_comp.dcenter
	direction = parse_direction(direction)
	if direction==1:# West edge
		distance = center.x - compbbox[0][0]
	elif direction==2:# North edge
		distance = center.y - compbbox[1][1]
	elif direction==3:# East edge
		distance = center.x - compbbox[1][0]
	elif direction==4:# South edge
		distance = center.y - compbbox[0][1]
	else:
		raise ValueError("unknown error with direction in function center_to_edge_distance (comp_utils)")
	return snap_to_grid(abs(distance),2)

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def move(custom_comp: Union[HPort, ComponentReference, Component], offsetxy: tuple[float,float] = (0,0), destination: Optional[tuple[Optional[float],Optional[float]]]=None, layer: Optional[tuple[int,int]]=None) -> Union[HPort, ComponentReference, Component]:
	"""moves custom_comp
	moves by offset[0]=x offset, offset[1]=y offset
	destination (x,y) if not none overrides offset option
	layer if specfied will move based on a layer (only relevant for destination option)
	returns the modified custom_comp
	"""
	if layer and isinstance(custom_comp, Component):
		custom_comp_ext = custom_comp.extract(layers=[layer,])
	elif layer and isinstance(custom_comp, ComponentReference):
		raise NotImplementedError("layer not implemented for comp ref")
	elif layer and isinstance(custom_comp,HPort):
		raise TypeError("move:layer option for HPort does not exist")
	else:
		custom_comp_ext = custom_comp
	if destination is not None:
		xoffset = destination[0] - (custom_comp_ext.dcenter.x) if destination[0] is not None else 0
		yoffset = destination[1] - (custom_comp_ext.dcenter.y) if destination[1] is not None else 0
	if isinstance(custom_comp, HPort):
		if destination is None:
			custom_comp = custom_comp.move_copy(offsetxy)
		else:
			custom_comp = custom_comp.move_copy((xoffset,yoffset))
	elif isinstance(custom_comp, ComponentReference):
		if destination is None:
			custom_comp.dmovex(offsetxy[0]).dmovey(offsetxy[1])
		else:
			custom_comp.dmovex(xoffset).dmovey(yoffset)
	elif isinstance(custom_comp, Component):
		ref = custom_comp.copy().add_ref(custom_comp)
		if destination is None:
			ref.dmovex(offsetxy[0]).dmovey(offsetxy[1])
		else:
			ref.dmovex(xoffset).dmovey(yoffset)
		custom_comp = transformed(ref).copy()
	return custom_comp


@validate_arguments(config=dict(arbitrary_types_allowed=True))
def movex(custom_comp: Union[HPort, ComponentReference, Component], offsetx: Optional[float] = 0, destination: Optional[float]=None, layer: Optional[tuple[int,int]]=None) -> Union[HPort, ComponentReference, Component]:
	"""moves custom_comp by offsetx in the x direction
	returns the modified custom_comp
	"""
	if destination is not None:
		destination = (destination, None)
	return move(custom_comp, (offsetx,0),destination,layer)


@validate_arguments(config=dict(arbitrary_types_allowed=True))
def movey(custom_comp: Union[HPort, ComponentReference, Component], offsety: Optional[float] = 0, destination: Optional[float]=None, layer: Optional[tuple[int,int]]=None) -> Union[HPort, ComponentReference, Component]:
	"""moves custom_comp by offsety in the y direction
	returns the modified custom_comp
	"""
	if destination is not None:
		destination = (None, destination)
	return move(custom_comp, (0,offsety), destination, layer)


# @validate_arguments(config=dict(arbitrary_types_allowed=True))
def align_comp_to_port(
	custom_comp: Union[Component,ComponentReference],
	align_to: HPort,
	alignment: Optional[tuple[Optional[str],Optional[str]]] = None,
	layer: Optional[tuple[int,int]] = None,
	rtr_comp_ref = True,
	port_width: Optional[float] = None,
	port_center: Optional[tuple[float,float]] = None, 
	return_only_coords: Optional[bool] = False
) -> Union[Component,ComponentReference]:
	"""Returns component/componentReference of component/componentReference aligned to port as specifed
	by default returns a componentReference
	for componentReference, the componentReference is modified (mutable), but for component, a copy of component is returned
	args:
	custom_comp = component to align properly
	align_to = HPort to align to
	alignment = tuple(str,str) = (xalign,yalign). You can individually specify x/y algin=None and that means do nothing for that dim
	***NOTE, if left None, function will align component to outside and center of port (based on port orientation), specify (None,None) for real no align (do not move at all)
	****xalign = either l/left or r/right or c/center or None. component will be flush to right or left side of port or centered
	****yalgin = either t/top or b/bottom or c/center or None. top or bottom edge or center of component will align with port top/bottom/center
	layer = extract this layer from the component and aligns to this layer.
	rtr_comp_ref = will return a component reference if set true, else return component
	"""
	# find center and bbox
	if isinstance(custom_comp, ComponentReference):
		comp_type = transformed(custom_comp)
	else:
		comp_type = custom_comp
	if layer:
		comp_type = comp_type.extract([layer])
	cbbox = comp_type.bbox_np()
	ccenter = comp_type.dcenter
	# setup
	# import pdb; pdb.set_trace()
	xdim = abs(cbbox[1][0] - cbbox[0][0])
	ydim = abs(cbbox[1][1] - cbbox[0][1])
	
	width = align_to.dwidth if not port_width else port_width
	alcenter = tuple(cc for cc in align_to.dcenter) if not port_center else port_width
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
	x_movcenter = alignto_center[0] - ccenter.x
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
	y_movcenter = alignto_center[1] - ccenter.y
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

	if return_only_coords: 
		retval = (xmov, ymov)
		return retval
	# make reference type, execute move
	if isinstance(custom_comp, Component):
		comp_ref = custom_comp.copy().add_ref(custom_comp)
	else:
		comp_ref = custom_comp
	# xmov = xmov * 1e-3
	# ymov = ymov * 1e-3
	# print(xmov, ymov)
	# import pdb; pdb.set_trace()
	comp_ref.dmovex(xmov).dmovey(ymov)
	# custom_comp = move(custom_comp, destination=(float(xmov),float(ymov)))
	# make correct type and return
	
	if rtr_comp_ref:
		return comp_ref
	else:
		return transformed(comp_ref)


@validate_arguments(config=dict(arbitrary_types_allowed=True))
def to_decimal(elements: Union[tuple,list,float,int,str]):
	"""converts all elements of list like object into decimals
	or converts single num into decimal"""
	if not isinstance(elements,Iterable):
		return Decimal(str(elements))
	else:
		elements = list(elements)
	for i, element in enumerate(elements):
		if isinstance(element,Union[int,float]):
			elements[i] = Decimal(str(element))
	return elements

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def to_float(elements: Union[tuple,list,Decimal,float]):
	"""converts all elements of list like object into floats and snaps to grid
	or converts single decimal into floats"""
	if not isinstance(elements,Iterable):
		return snap_to_grid(float(elements))
	else:
		elements = list(elements)
	for i, element in enumerate(elements):
		if isinstance(element, Union[float,Decimal]):
			elements[i] = snap_to_grid(float(element))
	return elements

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def prec_array(custom_comp: Component, rows: int, columns: int, spacing: tuple[Union[float,Decimal],Union[float,Decimal]], absolute_spacing: Optional[bool]=False) -> Component:
	"""instead of using the component.add_array function, if you are having grid snapping issues try using this function
	works the same way as add_array but uses decimals and snaps to grid to mitigate grid snapping issues
	args
	custom_comp: Component type to make an array from
	columns: num cols in the array
	rows: num rows in the array
	absolute_spacing: the spacing mode of spacing variable
	spacing: IF absolute_spacing spacing BETWEEN elements in the array ELSE spacing BETWEEN ORIGINS of elements in the array
	****NOTE do not use negative spacing, instead specify absolute_spacing=True
	"""
	# make sure to work with decimals
	precspacing = list(spacing)
	# import pdb; pdb.set_trace()
	
	for i in range(2):
		if isinstance(spacing[i],Union[int,float]):
			precspacing[i] = Decimal(str(spacing[i]))
	if not absolute_spacing:
		precspacing = [precspacing[i] + evaluate_bbox(custom_comp,True)[i] for i in range(2)]
	# create array
	precarray = Component()
	for colnum in range(columns):
		coldisp = colnum * precspacing[0]
		for rownum in range(rows):
			rowdisp = rownum * precspacing[1]
			cref = precarray << custom_comp
			cref.dmovex(to_float(coldisp)).dmovey(to_float(rowdisp))
			precarray.add_ports(cref.ports,prefix=f"row{rownum}_col{colnum}_")
	# precarray.flatten()
	return precarray

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def prec_center(custom_comp: Union[Component,ComponentReference], return_decimal: bool=False) -> tuple[Union[float,Decimal],Union[float,Decimal]]:
	"""instead of using component.ref_center() to get the center of a component,
	use this function which will return the correct offset to center a component
	returns (x,y) corrections
	if return_decimal=True, return in Decimal, otherwise return float"""
	
	correctmax = [dim/2 for dim in evaluate_bbox(custom_comp, True)]
	currentmax = to_decimal((custom_comp.dxmax,custom_comp.dymax))
	correctionxy = [correctmax[i] - currentmax[i] for i in range(2)]
	if return_decimal:
		return correctionxy
	return to_float(correctionxy)


@validate_arguments(config=dict(arbitrary_types_allowed=True))
def prec_ref_center(custom_comp: Union[Component,ComponentReference], destination: Optional[tuple[float,float]]=None, snapmov2grid: bool=False, return_comp: Optional[bool] = False):
	"""instead of using component.ref_center() to get a ref to center at origin,
	use this function which will return a centered ref
	you can then run component.add(prec_ref_center(custom_comp)) to add the reference to your component
	you can optionally specify a destination = tuple(x,y), by default, the destination is 0,0
	returns component reference
	"""
	xcor, ycor = prec_center(custom_comp, False)
	c = Component()
	compref = custom_comp if isinstance(custom_comp, ComponentReference) else c.add_ref(custom_comp)
	if destination is not None:
		xcor += destination[0]
		ycor += destination[1]
	if snapmov2grid:
		compref = movex(compref, snap_to_grid(xcor, 2))
		compref = movey(compref, snap_to_grid(ycor, 2))
	else:
		compref = movex(compref, xcor)
		compref = movey(compref, ycor)
		return compref
		# return compref.movex(xcor).movey(ycor)



def get_padding_points_cc(
	custom_comp: Union[ComponentReference, Component, list],
	default: float = 50.0,
	top: Optional[float]=None,
	bottom: Optional[float]=None,
	right: Optional[float]=None,
	left: Optional[float]=None,
	pdk_for_snap2xgrid: Optional[MappedPDK]=None
) -> list:
	"""works like gdsfactory.add_padding.get_padding_points, but also accepts componentReference or bbox
	additionally, if you optionally pass a pdk it will snap to 2x grid (else just operates like get_padding_points)"""
	if isinstance(custom_comp, Component):
		bbox = custom_comp.bbox_np()
	elif isinstance(custom_comp, ComponentReference):
		bbox = transformed(custom_comp).bbox_np()
	else:
		bbox = custom_comp
	top = top if top is not None else default
	bottom = bottom if bottom is not None else default
	right = right if right is not None else default
	left = left if left is not None else default
	ppoints = [
		[bbox[0][0] - left, bbox[0][1] - bottom],
		[bbox[1][0] + right, bbox[0][1] - bottom],
		[bbox[1][0] + right, bbox[1][1] + top],
		[bbox[0][0] - left, bbox[1][1] + top],
	]
	if pdk_for_snap2xgrid is not None:
		for i, ppoint in enumerate(ppoints):
			ppoints[i] = pdk_for_snap2xgrid.snap_to_2xgrid(ppoint)
	return ppoints

# function to take a component, delete duplicate ports
def delete_duplicate_ports(comp: Component) -> Component: 
    # get all ports
	ports = comp.ports
	# create a dictionary to store the ports
	port_dict = {}
	# iterate through all the ports
	for port in ports:
		# get the name of the port
		port_name = port.name
		# if the port name is not in the dictionary, add it
		if port_name not in port_dict:
			port_dict[port_name] = port
	# create a new component
	new_comp = Component()
	# add the ports to the new component
	new_comp.add_ports(port_dict.values())
	# return the new component
	return new_comp

def get_primitive_rectangle(size: tuple[float,float]=(5,3), dest: Optional[tuple[float, float]] = None, layer: tuple[int,int]=(0,0)):
	"""creates a rectangle component which snaps point to grid (does not snap to 2x grid)
	has same behavoir as gdsfactory.components.rectangle but doesnt allow centering (would snap to grid)
	"""
	c = Component()
	dx, dy = snap_to_grid2x(size)
	centered = False
	# import pdb; pdb.set_trace()
	if dx <= 0 or dy <= 0:
		raise ValueError(f"dx={dx} and dy={dy} must be > 0")
	
	if dest is None:
		points = [
			[0, 0],
			[0, dy],
			[dx, dy],
			[dx, 0],
		]
	else:
		points = [
			[dest[0], dest[1]],
			[dest[0], dest[1] + dy], 
			[dest[0] + dx, dest[1] + dy],
   			[dest[0] + dx, dest[1]],	
		]

	c.add_polygon(points, layer=layer)

	port_inclusion = 0.0
	port_type = "electrical"
	c.add_port(name="e1",center=(0, dy / 2),width=dy,orientation=180,layer=layer,port_type=port_type,)
	c.add_port(name="e2",center=(dx / 2 , dy),width=dx,orientation=90,layer=layer,port_type=port_type,)
	c.add_port(name="e3",center=(dx, dy / 2),width=dy,orientation=0,layer=layer,port_type=port_type,)
	c.add_port(name="e4",center=(dx / 2, 0),width=dx,orientation=-90,layer=layer,port_type=port_type,)
	return c
	# from gdsfactory.pdk import get_layer
	# temprect = Component()
	# # layer1 = get_layer(layer)
	# # print(f'\n\n{layer1}\n\n')
	# import pdb; pdb.set_trace()
	# temprect = rectangle(size, layer, centered = True, port_orientations=)
	# # shape = temprect.add_polygon(rectangle(size, layer, centered = True).bbox_np(), layer=layer)
	# # temprect.add_ports(temprect.get_ports_list(), prefix="route_")
	# temprect = add_ports_perimeter(temprect,layer=layer,prefix="route_")
	# temprect = rename_ports_by_list(temprect,(("W","e1"),("N","e2"),("E","e3"),("S","e4")))
	# # temprect.auto_rename_ports(temprect.ports,(("W","e1"),("N","e2"),("E","e3"),("S","e4")))
	# #rect = Component()
	# #clogic_ref = prec_ref_center(temprect) if centered else temprect.ref()
	# #rect.add(clogic_ref)
	# #rect.add_ports(clogic_ref.ports)
	# temprect = delete_duplicate_ports(temprect)
	# # temprect.flatten()
	# return temp


def create_rectangular_ring(
    enclosed_rectangle = (4.0, 2.0), 
    ring_width: float = 0.5, 
    layer: str = None, 
    centered: bool = False
) -> Component:
    if layer is None: 
        raise ValueError("layer must be specified!")
    if type(enclosed_rectangle) == list:
        # convert to tuple
        enclosed_rectangle = tuple(enclosed_rectangle)
    c = Component()
    inner = c.add_ref(rectangle(size=enclosed_rectangle, centered=centered, layer=layer))
    outerdims = (enclosed_rectangle[0] + 2 * ring_width, enclosed_rectangle[1] + 2 * ring_width)
    outer = rectangle(size=outerdims, centered=centered, layer=layer)
    if not centered:
        inner.move((ring_width, ring_width))
        
    return boolean(A = outer, B = inner, operation = 'not', layer=layer)

from gdsfactory import cell
@cell(check_ports = False)
def component_array(
    component: Component,
    spacing: tuple[float, float] = (150.0, 150.0),
    columns: int = 6,
    rows: int = 1,
    add_ports: bool = True,
    size: tuple[float, float] | None = None,
    centered: bool = False,
) -> Component:
    if size:
        columns = int(size[0] / spacing[0])
        rows = int(size[1] / spacing[1])

    if rows > 1 and spacing[1] == 0:
        raise ValueError(f"rows = {rows} > 1 require spacing[1] > 0")

    if columns > 1 and spacing[0] == 0:
        raise ValueError(f"columns = {columns} > 1 require spacing[0] > 0")

    c = Component()
    ref = c.add_ref(component, columns=columns, rows=rows, spacing=spacing)
    old_center = ref.dcenter
    ref.dcenter = (0, 0) if centered else old_center

    if add_ports and component.ports:
        for ix in range(ref.na):
            for iy in range(ref.nb):
                for port in component.ports:
                    port = port.copy(ref.trans * gf.kdb.Trans(ix * ref.a + iy * ref.b))
                    name = f"row{iy}_col{ix}_{port.name}"
                    c.add_port(name, port=port)

    return c