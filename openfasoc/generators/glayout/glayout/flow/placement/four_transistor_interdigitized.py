# 4 transistor placed in two rows (each row is an interdigitized pair of transistors)
# the 4 transistors are labeled top or bottom and transistor A or B
# top_A_, bottom_A, top_B_, bottom_B_

from glayout.flow.pdk.mappedpdk import MappedPDK
from glayout.flow.placement.two_transistor_interdigitized import two_nfet_interdigitized, two_pfet_interdigitized
from typing import Literal, Optional
from gdsfactory import Component
from glayout.flow.pdk.util.comp_utils import evaluate_bbox, movey
from glayout.flow.primitives.guardring import tapring
from gdsfactory.snap import snap_to_grid

def generic_4T_interdigitzed(
    pdk: MappedPDK,
    top_row_device: Literal["nfet", "pfet"],
    bottom_row_device: Literal["nfet", "pfet"],
    numcols: int,
    length: float=None,
    with_substrate_tap: bool = True,
    top_kwargs: Optional[dict]=None,
    bottom_kwargs: Optional[dict]=None
):
    if top_kwargs is None:
        top_kwargs = dict()
    if bottom_kwargs is None:
        bottom_kwargs = dict()
    # place
    toplvl = Component()
    if top_row_device=="nfet":
        toprow = toplvl << two_nfet_interdigitized(pdk,numcols,with_substrate_tap=False,length=length,**top_kwargs)
    else:
        toprow = toplvl << two_pfet_interdigitized(pdk,numcols,with_substrate_tap=False,length=length,**top_kwargs)
    toplvl.add_ports(toprow.ports,prefix="top_")
    
    if bottom_row_device=="nfet":
        bottomrow = toplvl << two_nfet_interdigitized(pdk,numcols,with_substrate_tap=False,length=length,**bottom_kwargs)
    else:
        bottomrow = toplvl << two_pfet_interdigitized(pdk,numcols,with_substrate_tap=False,length=length,**bottom_kwargs)
    toplvl.add_ports(bottomrow.ports,prefix="bottom_")

    # move
    toprow.dmovey(snap_to_grid((0.5 * (evaluate_bbox(bottomrow)[1]/2 + evaluate_bbox(toprow)[1]/2 + pdk.util_max_metal_seperation())), nm = 10))
    bottomrow.dmovey(-1 * snap_to_grid((0.5 * (evaluate_bbox(bottomrow)[1]/2 + evaluate_bbox(toprow)[1]/2 + pdk.util_max_metal_seperation())), nm = 10))
    # add substrate tap
    if with_substrate_tap:
        toplvl.flatten()
        substrate_tap_ref = toplvl << tapring(pdk, enclosed_rectangle=snap_to_grid(evaluate_bbox(toplvl,padding=pdk.util_max_metal_seperation()), nm = 10))
    # add ports
    toplvl.add_ports(substrate_tap_ref.ports,prefix="substratetap_")
    # flag for smart route
    toplvl.info["route_genid"] = "four_transistor_interdigitized"
    return toplvl
