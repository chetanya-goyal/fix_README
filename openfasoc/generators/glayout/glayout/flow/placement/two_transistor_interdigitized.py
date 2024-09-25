from glayout.flow.pdk.mappedpdk import MappedPDK
from pydantic import validate_arguments
from gdsfactory.component import Component
from glayout.flow.primitives.fet import nmos, pmos, multiplier
from glayout.flow.pdk.util.comp_utils import evaluate_bbox
from typing import Literal, Union
from glayout.flow.pdk.util.port_utils import rename_ports_by_orientation, rename_ports_by_list, create_private_ports
from glayout.flow.pdk.util.comp_utils import prec_ref_center, move, movex, movey
from glayout.flow.routing.straight_route import straight_route
from glayout.flow.pdk.util.comp_utils import transformed
from glayout.flow.primitives.guardring import tapring
from glayout.flow.pdk.util.port_utils import add_ports_perimeter
from gdsfactory import clear_cache
from typing import Literal
from gdsfactory.snap import snap_to_grid
from glayout.flow.routing.hashport import init_hashport
from gdsfactory.add_padding import add_padding

@validate_arguments(config=dict(arbitrary_types_allowed=True))
def macro_two_transistor_interdigitized(
    pdk: MappedPDK,
    numcols: int,
    deviceA_and_B: Literal["nfet", "pfet"],
    dummy: Union[bool, tuple[bool, bool]] = True,
    **kwargs
) -> Component:
    """place two transistors in a single row with interdigitized placement
    Currently only supports two of the same transistor (same devices)
    Place follows an ABABAB... pattern
    args:
    pdk = MappedPDK to use
    numcols = a single col is actually one col for both transistors (so AB). 2 cols = ABAB ... so on
    deviceA_and_B = the device to place for both transistors (either nfet or pfet)
    dummy = place dummy at the edges of the interdigitized place (true by default). you can specify tuple to place only on one side
    kwargs = key word arguments for device. 
    ****NOTE: These are the same as glayout.flow.primitives.fet.multiplier arguments EXCLUDING dummy, sd_route_extension, and pdk options
    """
    if isinstance(dummy, bool):
        dummy = (dummy, dummy)
    # override kwargs for needed options
    kwargs["sd_route_extension"] = 0
    kwargs["gate_route_extension"] = 0
    kwargs["sdlayer"] = "n+s/d" if deviceA_and_B == "nfet" else "p+s/d"
    kwargs["pdk"] = pdk
    # create devices dummy l/r and A/B (change extension options)
    kwargs["dummy"] = (True,False) if dummy[0] else False
    lefttmost_devA = multiplier(**kwargs)
    kwargs["dummy"] = False
    center_devA = multiplier(**kwargs)
    devB_sd_extension = pdk.util_max_metal_seperation() + abs(center_devA.ports["drain_N"].dcenter[1]-center_devA.ports["diff_N"].dcenter[1])
    devB_gate_extension = pdk.util_max_metal_seperation() + abs(center_devA.ports["row0_col0_gate_S"].dcenter[1]-center_devA.ports["gate_S"].dcenter[1])
    kwargs["sd_route_extension"] = snap_to_grid(0.8 * devB_sd_extension, nm = 10)
    kwargs["gate_route_extension"] = snap_to_grid(devB_gate_extension, nm = 10)
    center_devB = multiplier(**kwargs)
    kwargs["dummy"] = (False,True) if dummy[1] else False
    rightmost_devB = multiplier(**kwargs)
    
    # place devices
    idplace = Component()
    dims = evaluate_bbox(center_devA)
    xdisp = snap_to_grid(dims[0]+pdk.get_grule("active_diff")["min_separation"], nm = 10)
    refs = list()
    for i in range(2*numcols):
        if i==0:
            refs.append(idplace << lefttmost_devA)
        elif i==((2*numcols)-1):
            refs.append(idplace << rightmost_devB)
        elif i%2: # is odd (so device B)
            refs.append(idplace << center_devB)
        else: # not i%2 == i is even (so device A)
            refs.append(idplace << center_devA)
        refs[-1].dmovex(snap_to_grid(i*(xdisp), nm = 10))
        devletter = "B" if i%2 else "A"
        prefix=devletter+"_"+str(int(i/2))+"_"
        idplace.add_ports(refs[-1].ports, prefix=prefix)
    # extend poly layer for equal parasitics
    for i in range(2*numcols):
        desired_end_layer = 'met1'
        # desired_end_layer = pdk.layer_to_glayer(refs[i].ports["row0_col0_rightsd_top_met_N"].layer)
        idplace << straight_route(pdk, init_hashport(refs[i].ports["row0_col0_rightsd_top_met_N"], layer_overload=pdk.get_glayer('met1')),init_hashport(refs[-1].ports["drain_E"], layer_overload=pdk.get_glayer('met2')),glayer2=desired_end_layer, width = refs[0].ports["source_W"].dwidth * 0.8)
        idplace << straight_route(pdk, init_hashport(refs[i].ports["leftsd_top_met_N"], layer_overload=pdk.get_glayer('met1')),init_hashport(refs[-1].ports["drain_E"], layer_overload = pdk.get_glayer('met2')),glayer2=desired_end_layer, width = refs[0].ports["source_W"].dwidth * 0.8)
        if not i%2:
            desired_gate_end_layer = "poly"
            idplace << straight_route(pdk, init_hashport(refs[i].ports["row0_col0_gate_S"], layer_overload = pdk.get_glayer('poly')), init_hashport(refs[-1].ports["gate_E"], layer_overload = pdk.get_glayer('met2')),glayer2=desired_gate_end_layer)
    # merge s/d layer for all transistors
    idplace << straight_route(pdk, init_hashport(refs[0].ports["plusdoped_W"], layer_overload = pdk.get_glayer('n+s/d')),init_hashport(refs[-1].ports["plusdoped_E"], layer_overload = pdk.get_glayer('n+s/d')))
    # create s/d/gate connections extending over entire row
    A_src = idplace << rename_ports_by_orientation(rename_ports_by_list(straight_route(pdk, init_hashport(refs[0].ports["source_W"], layer_overload = pdk.get_glayer('met2')), init_hashport(refs[-1].ports["source_E"], layer_overload=pdk.get_glayer('met2')), width = refs[0].ports["source_W"].dwidth * 0.8), [("route_","_")]))
    B_src = idplace << rename_ports_by_orientation(rename_ports_by_list(straight_route(pdk, init_hashport(refs[-1].ports["source_E"], layer_overload=pdk.get_glayer('met2')), init_hashport(refs[0].ports["source_W"], layer_overload=pdk.get_glayer('met2')),width = refs[0].ports["source_W"].dwidth * 0.8), [("route_","_")]))
    A_drain = idplace << rename_ports_by_orientation(rename_ports_by_list(straight_route(pdk, init_hashport(refs[0].ports["drain_W"], layer_overload=pdk.get_glayer('met2')), init_hashport(refs[-1].ports["drain_E"], layer_overload=pdk.get_glayer('met2')),width = refs[0].ports["source_W"].dwidth * 0.8), [("route_","_")]))
    B_drain = idplace << rename_ports_by_orientation(rename_ports_by_list(straight_route(pdk, init_hashport(refs[-1].ports["drain_E"], layer_overload=pdk.get_glayer('met2')), init_hashport(refs[0].ports["drain_W"], layer_overload=pdk.get_glayer('met2')),width = refs[0].ports["source_W"].dwidth * 0.8), [("route_","_")]))
    A_gate = idplace << rename_ports_by_orientation(rename_ports_by_list(straight_route(pdk, init_hashport(refs[0].ports["gate_W"], layer_overload=pdk.get_glayer('met2')), init_hashport(refs[-1].ports["gate_E"], layer_overload=pdk.get_glayer('met2'))), [("route_","_")]))
    B_gate = idplace << rename_ports_by_orientation(rename_ports_by_list(straight_route(pdk, init_hashport(refs[-1].ports["gate_E"], layer_overload=pdk.get_glayer('met2')), init_hashport(refs[0].ports["gate_W"], layer_overload=pdk.get_glayer('met2'))), [("route_","_")]))
    # add route ports and return

    prefixes = ["A_source","B_source","A_drain","B_drain","A_gate","B_gate"]
    for i, ref in enumerate([A_src, B_src, A_drain, B_drain, A_gate, B_gate]):
        idplace.add_ports(ref.ports,prefix=prefixes[i])
    idplace = transformed(prec_ref_center(idplace))
    # idplace.unlock()
    idplace.add_ports(create_private_ports(idplace, prefixes))
    return idplace


@validate_arguments(config=dict(arbitrary_types_allowed=True))
def two_nfet_interdigitized(
    pdk: MappedPDK,
    numcols: int,
    dummy: Union[bool, tuple[bool, bool]] = True,
    with_substrate_tap: bool = True,
    with_tie: bool = True,
    tie_layers: tuple[str,str]=("met2","met1"),
    **kwargs
) -> Component:
    """Currently only supports two of the same nfet instances. does NOT support multipliers (currently)
    Place follows an ABABAB... pattern
    args:
    pdk = MappedPDK to use
    numcols = a single col is actually one col for both nfets (so AB). 2 cols = ABAB ... so on
    dummy = place dummy at the edges of the interdigitized place (true by default). you can specify tuple to place only on one side
    kwargs = key word arguments for multiplier. 
    ****NOTE: These are the same as glayout.flow.primitives.fet.multiplier arguments EXCLUDING dummy, sd_route_extension, and pdk options
    tie_layers: tuple[str,str] specifying (horizontal glayer, vertical glayer) or well tie ring. default=("met2","met1")
    """
    base_multiplier = Component()
    dummycomp = base_multiplier << macro_two_transistor_interdigitized(pdk, numcols, "nfet", dummy, **kwargs)
    base_multiplier.add_ports(dummycomp.ports)
    dummycomp.dmove(origin = (dummycomp.dcenter.x, dummycomp.dcenter.y), destination = (0, 0))
    # tie
    
    if with_tie:
        tap_separation = max(
            pdk.util_max_metal_seperation(),
            pdk.get_grule("active_diff", "active_tap")["min_separation"],
        )
        tap_separation += pdk.get_grule("p+s/d", "active_tap")["min_enclosure"]
        tap_encloses = (
            snap_to_grid(2 * (tap_separation + base_multiplier.dxmax), nm = 10),
            snap_to_grid(2 * (tap_separation + base_multiplier.dymax), nm = 10)
        )
        tiering_ref = base_multiplier << tapring(
            pdk,
            enclosed_rectangle=tap_encloses,
            sdlayer="p+s/d",
            horizontal_glayer=tie_layers[0],
            vertical_glayer=tie_layers[1],
        )
        base_multiplier.add_ports(tiering_ref.ports, prefix="welltie_")
        # try:
        #     base_multiplier<<straight_route(pdk,init_hashport(base_multiplier.ports["A_0_dummy_L_gsdcon_top_met_W"], layer_overload = pdk.get_glayer('met1')), init_hashport(base_multiplier.ports["welltie_W_top_met_W"], layer_overload = pdk.get_glayer('met1')), glayer2="met1")
        # except KeyError:
        #     pass
        # try:
        #     base_multiplier<<straight_route(pdk,init_hashport(base_multiplier.ports[f"B_{numcols-1}_dummy_R_gsdcon_top_met_E"], layer_overload = pdk.get_glayer('met1')), init_hashport(base_multiplier.ports["welltie_E_top_met_E"], layer_overload = pdk.get_glayer('met1')), glayer2="met1")
        # except KeyError:
        #     pass
    # add pwell
    base_multiplier = add_padding(
        base_multiplier,
        layers=(pdk.get_glayer("pwell"),),
        default=pdk.get_grule("pwell", "active_tap")["min_enclosure"],
    )
    # add substrate tap
    base_multiplier = add_ports_perimeter(base_multiplier,layer=pdk.get_glayer("pwell"),prefix="well_")
    # add substrate tap if with_substrate_tap
    if with_substrate_tap:
        substrate_tap_separation = pdk.get_grule("dnwell", "active_tap")[
            "min_separation"
        ]
        substrate_tap_encloses = (
            snap_to_grid(2 * (substrate_tap_separation + base_multiplier.dxmax), nm = 10),
            snap_to_grid(2 * (substrate_tap_separation + base_multiplier.dymax), nm = 10)
        )
        ringtoadd = tapring(
            pdk,
            enclosed_rectangle=substrate_tap_encloses,
            sdlayer="p+s/d",
            horizontal_glayer="met2",
            vertical_glayer="met1",
        )
        tapring_ref = base_multiplier << ringtoadd
        base_multiplier.add_ports(tapring_ref.ports,prefix="substratetap_")
    base_multiplier.info["route_genid"] = "two_transistor_interdigitized"
    return base_multiplier



@validate_arguments(config=dict(arbitrary_types_allowed=True))
def two_pfet_interdigitized(
    pdk: MappedPDK,
    numcols: int,
    dummy: Union[bool, tuple[bool, bool]] = True,
    with_substrate_tap: bool = True,
    with_tie: bool = True,
    tie_layers: tuple[str,str]=("met2","met1"),
    **kwargs
) -> Component:
    """Currently only supports two of the same nfet instances. does NOT support multipliers (currently)
    Place follows an ABABAB... pattern
    args:
    pdk = MappedPDK to use
    numcols = a single col is actually one col for both nfets (so AB). 2 cols = ABAB ... so on
    dummy = place dummy at the edges of the interdigitized place (true by default). you can specify tuple to place only on one side
    kwargs = key word arguments for multiplier. 
    ****NOTE: These are the same as glayout.flow.primitives.fet.multiplier arguments EXCLUDING dummy, sd_route_extension, and pdk options
    tie_layers: tuple[str,str] specifying (horizontal glayer, vertical glayer) or well tie ring. default=("met2","met1")
    """
    base_multiplier = Component()
    dummycomp = base_multiplier << macro_two_transistor_interdigitized(pdk, numcols, "pfet", dummy, **kwargs)
    base_multiplier.add_ports(dummycomp.ports)
    dummycomp.dmove(origin = (dummycomp.dcenter.x, dummycomp.dcenter.y), destination = (0, 0))
    # tie
    if with_tie:
        tap_separation = max(
            pdk.util_max_metal_seperation(),
            pdk.get_grule("active_diff", "active_tap")["min_separation"],
        )
        tap_separation += pdk.get_grule("n+s/d", "active_tap")["min_enclosure"]
        tap_encloses = (
            snap_to_grid(2 * (tap_separation + base_multiplier.dxmax), nm = 10),
            snap_to_grid(2 * (tap_separation + base_multiplier.dymax), nm = 10)
        )
        tiering_ref = base_multiplier << tapring(
            pdk,
            enclosed_rectangle=tap_encloses,
            sdlayer="n+s/d",
            horizontal_glayer=tie_layers[0],
            vertical_glayer=tie_layers[1],
        )
        base_multiplier.add_ports(tiering_ref.ports, prefix="welltie_")
        # try:
        #     base_multiplier<<straight_route(pdk,init_hashport(base_multiplier.ports["A_0_dummy_L_gsdcon_top_met_W"], layer_overload = pdk.get_glayer('met1')), init_hashport(base_multiplier.ports["welltie_W_top_met_W"], layer_overload = pdk.get_glayer('met1')), glayer2="met1")
        # except KeyError:
        #     pass
        # try:
        #     base_multiplier<<straight_route(pdk, init_hashport(base_multiplier.ports[f"B_{numcols-1}_dummy_R_gsdcon_top_met_E"], layer_overload = pdk.get_glayer('met1')), init_hashport(base_multiplier.ports["welltie_E_top_met_E"], layer_overload = pdk.get_glayer('met1')), glayer2="met1")
        # except KeyError:
        #     pass
    # add pwell
    base_multiplier = add_padding(
        base_multiplier,
        layers=(pdk.get_glayer("nwell"),),
        default=pdk.get_grule("nwell", "active_tap")["min_enclosure"],
    )
    # add substrate tap
    base_multiplier = add_ports_perimeter(base_multiplier,layer=pdk.get_glayer("nwell"),prefix="well_")
    # add substrate tap if with_substrate_tap
    if with_substrate_tap:
        substrate_tap_separation = pdk.get_grule("dnwell", "active_tap")[
            "min_separation"
        ]
        substrate_tap_encloses = (
            snap_to_grid(2 * (substrate_tap_separation + base_multiplier.dxmax), nm = 10),
            snap_to_grid(2 * (substrate_tap_separation + base_multiplier.dymax), nm = 10)
        )
        ringtoadd = tapring(
            pdk,
            enclosed_rectangle=substrate_tap_encloses,
            sdlayer="p+s/d",
            horizontal_glayer="met2",
            vertical_glayer="met1",
        )
        tapring_ref = base_multiplier << ringtoadd
        base_multiplier.add_ports(tapring_ref.ports,prefix="substratetap_")
    base_multiplier.info["route_genid"] = "two_transistor_interdigitized"
    return base_multiplier




def two_transistor_interdigitized(
    pdk: MappedPDK,
    device: Literal["nfet","pfet"],
    numcols: int,
    dummy: Union[bool, tuple[bool, bool]] = True,
    with_substrate_tap: bool = True,
    with_tie: bool = True,
    tie_layers: tuple[str,str]=("met2","met1"),
    **kwargs
) -> Component:
    if device=="nfet":
        return two_nfet_interdigitized(pdk=pdk,numcols=numcols,dummy=dummy,with_substrate_tap=with_substrate_tap,with_tie=with_tie,tie_layers=tie_layers,**kwargs)
    else:
        return two_pfet_interdigitized(pdk=pdk,numcols=numcols,dummy=dummy,with_substrate_tap=with_substrate_tap,with_tie=with_tie,tie_layers=tie_layers,**kwargs)


