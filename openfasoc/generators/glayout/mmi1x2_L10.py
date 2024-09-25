from __future__ import annotations

import gdsfactory as gf

if __name__ == "__main__":
    c = (
        gf.Component()
    )
    p1 = c.add_polygon([(-8, -6), (6, 8), (7, 17), (9, 5)], layer=(1, 0))
    # c.write('demo.gds')
    c.write_gds('demo.gds')
    # c.show()
    # print(c.layers.index((1,0)))