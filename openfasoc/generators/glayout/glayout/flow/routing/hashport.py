from __future__ import annotations

import csv
import functools
import typing
import warnings
from collections.abc import Callable
from functools import partial
from typing import Any, Dict, List, Optional, Tuple, Union
import kfactory as kf
import numpy as np
from rich.console import Console
from rich.table import Table

from gdsfactory.cross_section import CrossSectionSpec
Layer = tuple[int, int]
Layers = tuple[Layer, ...]
LayerSpec = Layer | str | None | kf.LayerEnum
LayerSpecs = tuple[LayerSpec, ...]
Float2 = tuple[float, float]
valid_error_types = ["error", "warn", "ignore"]

class HPort(kf.Port):
    """Ports are useful to connect Components with each other.

    Args:
        name: we name ports clock-wise starting from bottom left.
        orientation: in degrees (0: east, 90: north, 180: west, 270: south).
        center: (x, y) port center coordinate.
        width: of the port in um.
        layer: layer tuple.
        port_type: str (optical, electrical, vertical_te, vertical_tm).
        cross_section: cross_section spec.
        info: additional information.
    """

    def __init__(
        self,
        name: str,
        orientation: float | None,
        center: tuple[float, float] | kf.kdb.Point | kf.kdb.DPoint,
        width: int,
        layer: LayerSpec | None = None,
        port_type: str = "optical",
        cross_section: CrossSectionSpec | None = None,
        info: dict[str, int | float | str] | None = None,
    ) -> None:
        """Initializes Port."""
        from gdsfactory.pdk import get_layer
        orientation = np.mod(orientation, 360) if orientation else orientation

        if cross_section is None and layer is None:
            raise ValueError("You need to define Port cross_section or layer")

        if cross_section is None and width is None:
            raise ValueError("You need Port to define cross_section or width")

        if layer is None or width is None:
            from gdsfactory.pdk import get_cross_section

            cross_section = get_cross_section(cross_section)

        if cross_section and layer is None:
            layer = cross_section.layer
        
        if isinstance(layer, list):
            layer = tuple(layer)

        if width is None:
            width = cross_section.width

        if width < 0:
            raise ValueError(f"Port width must be >=0. Got {width}")

        dcplx_trans = kf.kdb.DCplxTrans(kf.kcl.dbu, float(orientation), False, *center)
        info = info or {}
        width = width
        layer = layer if isinstance(layer, tuple) else get_layer(layer)
        super().__init__(
            name=name,
            layer=layer,
            width=width,
            port_type=port_type,
            trans=dcplx_trans.s_trans().to_itype(1.0)
        )
        
    def __eq__(self, other):
        """Check equality based on unique attributes."""
        if not isinstance(other, HPort):
            return NotImplemented
        return (self.name == other.name and self.orientation == other.orientation and
                self.center == other.center and self.width == other.width and
                self.layer == other.layer)

    def __hash__(self):
        """Generate hash based on unique attributes."""
        return hash((self.name, self.orientation, self.center, self.width, self.layer))
    
    def Hcopy(
        self, 
        name_overload: Optional[str] = None,
        orientation_overload: Optional[float] = None,
        center_overload: Optional[Union[tuple[Optional[float], Optional[float]], kf.kdb.Point, kf.kdb.DPoint]] = None,
        width_overload: Optional[int] = None,
        layer_overload: Optional[LayerSpec] = None
    ):
        """Creates a copy of the port with overloaded parameters."""
        name = self.name if name_overload is None else name_overload
        orientation = self.orientation if orientation_overload is None else orientation_overload
        center = self.center if center_overload is None else center_overload
        width = self.width if width_overload is None else width_overload
        layer = self.layer if layer_overload is None else layer_overload
        # newcenter = (int(float(center[0]) * kf.kcl.dbu), int(float(center[1]) * kf.kcl.dbu))
        return init_hashport(
            self, 
            name_overload=name,
            orientation_overload=orientation,
            center_overload=center,
            width_overload=width,
            layer_overload=layer
        )
        
    def move_copy(self, x: int = 0, y: int = 0):
        _port = self.Hcopy()
        _port.dcenter = (x, y)
        return _port
    
    def Hmovex(self, x: Union[int, float]):
        _port = self.Hcopy()
        _port.dcenter = (x, _port.dcenter[1])
        return _port
    
    def Hmovey(self, y: Union[int, float]):
        _port = self.Hcopy()
        _port.dcenter = (_port.dcenter[0], y)
        return _port

def init_hashport(
    input_port: kf.Port, 
    name_overload: Optional[str] = None,
    orientation_overload: Optional[float] = None,
    center_overload: Optional[Union[tuple[Optional[float], Optional[float]], kf.kdb.Point, kf.kdb.DPoint]] = None,
    width_overload: Optional[int] = None,
    layer_overload: Optional[LayerSpec] = None
) -> HPort:
    """creates a hashable port object from a kfactory port object"""
    name = input_port.name if name_overload is None else name_overload
    orientation = input_port.orientation if orientation_overload is None else orientation_overload
    center = input_port.center if center_overload is None else center_overload
    width = input_port.width if width_overload is None else width_overload
    layer = input_port.layer if layer_overload is None else layer_overload
    return HPort(
        name=name,
        orientation=orientation,
        center=center,
        width=width,
        layer=layer,
        port_type=input_port.port_type
    )