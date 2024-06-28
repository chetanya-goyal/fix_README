from gdsfactory.typings import Component
from pydantic import validate_arguments
from typing import Optional

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def component_snap_to_grid(comp: Component, flatten_true: Optional[bool] = True) -> Component:
	"""snaps all polygons and ports in component to grid
	comp = the component to snap to grid
	NOTE this function will flatten the component
	"""
	#return comp.flatten()
	# flatten the component then copy (snaps polygons and ports to grid)
	name = comp.name
	
	if flatten_true:
		comp.flatten()
	comp = comp.copy()
	comp.name = name
	return comp


