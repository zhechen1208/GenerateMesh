#!/usr/bin/env python3
"""
Convert Tecplot structured DAT files (POINT format) to a linear-hexahedral
Nektar++ XML mesh.

This is an initial/debug-friendly converter for structured O/H grids:
  - reads one or more Tecplot DAT files with `zone i=..., j=..., k=...`;
  - reconstructs implicit structured topology;
  - creates unique VERTEX, EDGE, FACE and HEXA elements;
  - handles periodic/closed j-seams automatically when j=1 and j=jmax coincide;
  - writes basic boundary COMPOSITE groups and a fluid DOMAIN;
  - optionally writes a simple 3D EXPANSIONS section.

Assumed point order in DAT:
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          write x, y, z

Typical use:
    python dat_to_nektar3d.py 3d_mesh.dat -o wing3d.xml
    python dat_to_nektar3d.py 3d_mesh.dat left_tip_block.dat right_tip_block.dat -o wing3d.xml --tol 1e-10

Notes:
  1. This writes linear hex elements: every structured cell becomes one Nektar++ <H>.
  2. For Nektar++ production calculations, a coarser macro-element mesh with curved
     faces is usually preferable.
  3. Vertex ordering uses j-flipped convention for right-handed hex elements.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

Point = Tuple[float, float, float]
EdgeKey = Tuple[int, int]
FaceKey = Tuple[int, int, int, int]


@dataclass
class Zone:
    name: str
    ni: int
    nj: int
    nk: int
    points: List[Point]
    source: str
    closed_j: bool = False

    def local_index(self, i: int, j: int, k: int) -> int:
        return k * self.ni * self.nj + j * self.ni + i

    def p(self, i: int, j: int, k: int) -> Point:
        return self.points[self.local_index(i, j, k)]


@dataclass
class Mesh3D:
    tol: float = 1.0e-10

    vertices: List[Point] = field(default_factory=list)
    edges: List[Tuple[int, int]] = field(default_factory=list)
    faces: List[Tuple[int, int, int, int]] = field(default_factory=list)  # edge IDs in loop order
    hexes: List[Tuple[int, int, int, int, int, int]] = field(default_factory=list)  # face IDs

    vertex_map: Dict[Tuple[int, int, int], int] = field(default_factory=dict)
    edge_map: Dict[EdgeKey, int] = field(default_factory=dict)
    face_map: Dict[FaceKey, int] = field(default_factory=dict)
    face_use_count: Dict[int, int] = field(default_factory=dict)
    face_tags: Dict[int, List[str]] = field(default_factory=dict)

    def coord_key(self, p: Point) -> Tuple[int, int, int]:
        if self.tol <= 0.0:
            # Exact key using raw floats is not hash-stable for intended merging;
            # still useful for disabling tolerance-based merge for debug.
            return (hash(p[0]), hash(p[1]), hash(p[2]))
        return (round(p[0] / self.tol), round(p[1] / self.tol), round(p[2] / self.tol))

    def add_vertex(self, p: Point) -> int:
        key = self.coord_key(p)
        if key in self.vertex_map:
            return self.vertex_map[key]
        vid = len(self.vertices)
        self.vertices.append(p)
        self.vertex_map[key] = vid
        return vid

    def add_edge(self, v0: int, v1: int) -> int:
        if v0 == v1:
            raise ValueError(f"Degenerate edge with identical vertices: {v0}")
        key = (v0, v1) if v0 < v1 else (v1, v0)
        eid = self.edge_map.get(key)
        if eid is not None:
            return eid
        eid = len(self.edges)
        self.edges.append((v0, v1))
        self.edge_map[key] = eid
        return eid

    def add_face_from_vertices(self, verts: Sequence[int], tag: Optional[str] = None) -> int:
        if len(verts) != 4:
            raise ValueError("Only quadrilateral faces are supported.")
        if len(set(verts)) != 4:
            raise ValueError(f"Degenerate face with repeated vertices: {verts}")

        # Face uniqueness is orientation-independent.
        key = tuple(sorted(verts))
        fid = self.face_map.get(key)
        if fid is None:
            e0 = self.add_edge(verts[0], verts[1])
            e1 = self.add_edge(verts[1], verts[2])
            e2 = self.add_edge(verts[2], verts[3])
            e3 = self.add_edge(verts[3], verts[0])
            fid = len(self.faces)
            self.faces.append((e0, e1, e2, e3))
            self.face_map[key] = fid
            self.face_use_count[fid] = 0
            self.face_tags[fid] = []
        self.face_use_count[fid] += 1
        if tag is not None:
            self.face_tags[fid].append(tag)
        return fid

    def add_hex(self, face_ids: Sequence[int]) -> int:
        if len(face_ids) != 6:
            raise ValueError("A hex must have six face IDs.")
        hid = len(self.hexes)
        self.hexes.append(tuple(face_ids))
        return hid


def parse_zone_header(line: str, default_name: str) -> Optional[Tuple[str, int, int, int]]:
    s = line.strip()
    if not s.lower().startswith("zone"):
        return None

    def grab_dim(name: str) -> Optional[int]:
        m = re.search(rf"\b{name}\s*=\s*(\d+)", s, flags=re.IGNORECASE)
        return int(m.group(1)) if m else None

    ni = grab_dim("i")
    nj = grab_dim("j")
    nk = grab_dim("k")
    if ni is None or nj is None:
        raise ValueError(f"Cannot parse zone dimensions from: {line}")
    if nk is None:
        nk = 1

    m = re.search(r't\s*=\s*"([^"]+)"', s, flags=re.IGNORECASE)
    name = m.group(1) if m else default_name
    return name, ni, nj, nk


def numeric_values(line: str) -> List[float]:
    out: List[float] = []
    for token in line.replace(",", " ").split():
        try:
            out.append(float(token))
        except ValueError:
            pass
    return out


def read_tecplot_dat(path: str) -> List[Zone]:
    p = Path(path)
    lines = p.read_text(errors="ignore").splitlines()
    zones: List[Zone] = []
    idx = 0
    zone_counter = 0

    while idx < len(lines):
        parsed = parse_zone_header(lines[idx], f"zone{zone_counter}")
        if parsed is None:
            idx += 1
            continue

        name, ni, nj, nk = parsed
        npts = ni * nj * nk
        idx += 1
        pts: List[Point] = []

        while idx < len(lines) and len(pts) < npts:
            vals = numeric_values(lines[idx])
            if len(vals) >= 3:
                pts.append((vals[0], vals[1], vals[2]))
            idx += 1

        if len(pts) != npts:
            raise RuntimeError(
                f"{path}: zone {name} expects {npts} points, only read {len(pts)}."
            )

        zones.append(Zone(name=name, ni=ni, nj=nj, nk=nk, points=pts, source=str(p)))
        zone_counter += 1

    if not zones:
        raise RuntimeError(f"No Tecplot structured zones found in {path}")
    return zones


def dist(a: Point, b: Point) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def detect_closed_j(zone: Zone, tol: float) -> bool:
    if zone.nj < 3:
        return False
    maxdiff = 0.0
    for k in range(zone.nk):
        for i in range(zone.ni):
            maxdiff = max(maxdiff, dist(zone.p(i, 0, k), zone.p(i, zone.nj - 1, k)))
    # Use a slightly relaxed criterion because DAT may have limited precision.
    return maxdiff <= max(10.0 * tol, 1.0e-12)


def compress_ids(ids: Iterable[int]) -> str:
    data = sorted(set(ids))
    if not data:
        return ""
    chunks: List[str] = []
    start = prev = data[0]
    for x in data[1:]:
        if x == prev + 1:
            prev = x
            continue
        chunks.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = x
    chunks.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(chunks)


def zone_is_tip_block(zone: Zone) -> bool:
    """Heuristic: separate standalone tip/cap block files from the main block."""
    s = (Path(zone.source).stem + "_" + zone.name).lower()
    return ("tip" in s) or ("cap" in s) or ("block" in s and ("left" in s or "right" in s))


def zone_i_min_tag(zone: Zone, k_cell: int, kleft: int | None, kright: int | None) -> str:
    """Return boundary tag for the i=1 face of a cell.

    k_cell is 0-based and denotes the interval between Tecplot node
    k=k_cell+1 and k=k_cell+2 in Fortran-style 1-based indexing.

    For the main block, when kleft/kright are provided, split the i=1
    wall into three stream/spanwise portions:
        cells [1, kleft]     -> main_left_transition_wall
        cells [kleft, kright]-> wall
        cells [kright, kmax] -> main_right_transition_wall
    More precisely, these are cell intervals:
        k_cell = 0 .. kleft-2       left transition
        k_cell = kleft-1 .. kright-2 wall
        k_cell = kright-1 .. nk-2    right transition
    """
    s = (Path(zone.source).stem + "_" + zone.name).lower()

    # Standalone tip/cap block files keep their own tags.
    if "left" in s and ("tip" in s or "block" in s):
        return "left_tip_wall"
    if "right" in s and ("tip" in s or "block" in s):
        return "right_tip_wall"
    if "tip" in s or "cap" in s:
        return "tip_wall"

    # Main block: split i=1 according to kleft/kright if given.
    if kleft is not None and kright is not None:
        if not (1 <= kleft < kright <= zone.nk):
            raise ValueError(
                f"Invalid kleft/kright for zone {zone.name}: "
                f"kleft={kleft}, kright={kright}, zone.nk={zone.nk}"
            )
        if k_cell <= kleft - 2:
            return "main_left_transition_wall"
        if k_cell >= kright - 1:
            return "main_right_transition_wall"
        return "wall"

    return "wall"


def build_zone(
    mesh: Mesh3D,
    zone: Zone,
    closed_j_auto: bool = True,
    force_open_j: bool = False,
    kleft: int | None = None,
    kright: int | None = None,
) -> None:
    zone.closed_j = False if force_open_j else (detect_closed_j(zone, mesh.tol) if closed_j_auto else False)

    print(
        f"Zone {zone.name} from {Path(zone.source).name}: "
        f"i={zone.ni}, j={zone.nj}, k={zone.nk}, closed_j={zone.closed_j}"
    )

    # Create local structured index -> global vertex ID.
    gids: Dict[Tuple[int, int, int], int] = {}
    j_vertex_count = zone.nj - 1 if zone.closed_j else zone.nj

    for k in range(zone.nk):
        for j in range(j_vertex_count):
            for i in range(zone.ni):
                gids[(i, j, k)] = mesh.add_vertex(zone.p(i, j, k))

    def gid(i: int, j: int, k: int) -> int:
        if zone.closed_j:
            j = j % (zone.nj - 1)
        return gids[(i, j, k)]

    # Cell range. For closed j, j=nj-2 connects to duplicate seam j=nj-1 -> j=0.
    j_cells = zone.nj - 1

    for k in range(zone.nk - 1):
        for j in range(j_cells):
            jp = j + 1
            if zone.closed_j and jp == zone.nj - 1:
                jp = 0
            for i in range(zone.ni - 1):
                # Vertex order with j-flipped to match Nektar++ right-handed convention.
                v000 = gid(i, j, k)
                v010 = gid(i, jp, k)
                v110 = gid(i + 1, jp, k)
                v100 = gid(i + 1, j, k)
                v001 = gid(i, j, k + 1)
                v011 = gid(i, jp, k + 1)
                v111 = gid(i + 1, jp, k + 1)
                v101 = gid(i + 1, j, k + 1)

                # Boundary tags are only used if the face is ultimately external.
                tag_kmin = "kmin" if k == 0 else None
                tag_kmax = "kmax" if k == zone.nk - 2 else None
                tag_jmin = "jmin" if (not zone.closed_j and j == 0) else None
                tag_jmax = "jmax" if (not zone.closed_j and j == zone.nj - 2) else None
                tag_imin = zone_i_min_tag(zone, k, kleft, kright) if i == 0 else None
                tag_imax = "farfield" if i == zone.ni - 2 else None

                # Six quad faces following Nektar++/Gmsh hex face convention.
                # Face 0: bottom (kmin), nodes 0-1-2-3
                f_kmin = mesh.add_face_from_vertices([v000, v100, v110, v010], tag_kmin)
                # Face 1: front (jmin), nodes 0-4-5-1
                f_jmin = mesh.add_face_from_vertices([v000, v001, v101, v100], tag_jmin)
                # Face 2: right (imax), nodes 1-5-6-2
                f_imax = mesh.add_face_from_vertices([v100, v101, v111, v110], tag_imax)
                # Face 3: back  (jmax), nodes 2-6-7-3
                f_jmax = mesh.add_face_from_vertices([v110, v111, v011, v010], tag_jmax)
                # Face 4: left  (imin), nodes 3-7-4-0
                f_imin = mesh.add_face_from_vertices([v010, v011, v001, v000], tag_imin)
                # Face 5: top   (kmax), nodes 4-7-6-5
                f_kmax = mesh.add_face_from_vertices([v001, v011, v111, v101], tag_kmax)

                # Nektar++ hex face order: bottom, front, right, back, left, top
                mesh.add_hex([f_kmin, f_jmin, f_imax, f_jmax, f_imin, f_kmax])


def boundary_groups(mesh: Mesh3D) -> Dict[str, List[int]]:
    """Collect external faces by their physical tag.

    A face is external if it is used by only one hex.  Tags are attached while
    building structured zones.  Unlike the first draft, this version keeps
    arbitrary tags such as left_tip_wall/right_tip_wall instead of putting them
    into "other".
    """
    preferred_order = [
        "wall",
        "main_left_transition_wall",
        "main_right_transition_wall",
        "left_tip_wall",
        "right_tip_wall",
        "tip_wall",
        "farfield",
        "kmin",
        "kmax",
        "jmin",
        "jmax",
        "other",
    ]

    groups: Dict[str, List[int]] = {name: [] for name in preferred_order}

    for fid in range(len(mesh.faces)):
        if mesh.face_use_count.get(fid, 0) != 1:
            continue

        tags = [t for t in mesh.face_tags.get(fid, []) if t is not None]
        if tags:
            # Use the first non-null tag.  For normal structured blocks a true
            # external face receives only one physical tag.
            name = tags[0]
        else:
            name = "other"

        if name not in groups:
            groups[name] = []
        groups[name].append(fid)

    ordered: Dict[str, List[int]] = {}
    for name in preferred_order:
        if groups.get(name):
            ordered[name] = groups[name]
    for name in sorted(groups.keys()):
        if name not in ordered and groups[name]:
            ordered[name] = groups[name]
    return ordered

def write_xml(
    path: str,
    mesh: Mesh3D,
    num_modes: int = 2,
    num_points: int = 3,
    fields: str = "u,v,w,p",
    write_expansions: bool = True,
) -> None:
    bgroups = boundary_groups(mesh)
    fluid_comp_id = len(bgroups)

    with open(path, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0" encoding="utf-8" ?>\n')
        f.write('<NEKTAR>\n')
        f.write('    <GEOMETRY DIM="3" SPACE="3">\n')

        f.write('        <VERTEX>\n')
        for i, (x, y, z) in enumerate(mesh.vertices):
            f.write(f'            <V ID="{i}">{x:26.17e}{y:26.17e}{z:26.17e}</V>\n')
        f.write('        </VERTEX>\n')

        f.write('        <EDGE>\n')
        for i, (v0, v1) in enumerate(mesh.edges):
            f.write(f'            <E ID="{i}">{v0} {v1} </E>\n')
        f.write('        </EDGE>\n')

        f.write('        <FACE>\n')
        for i, (e0, e1, e2, e3) in enumerate(mesh.faces):
            f.write(f'            <Q ID="{i}">{e0} {e1} {e2} {e3} </Q>\n')
        f.write('        </FACE>\n')

        f.write('        <ELEMENT>\n')
        for i, h in enumerate(mesh.hexes):
            f.write(
                f'            <H ID="{i}">{h[0]} {h[1]} {h[2]} '
                f'{h[3]} {h[4]} {h[5]} </H>\n'
            )
        f.write('        </ELEMENT>\n')

        f.write('        <CURVED />\n')

        f.write('        <COMPOSITE>\n')
        comp_id = 0
        for name, ids in bgroups.items():
            f.write(f'            <!-- {name} -->\n')
            f.write(f'            <C ID="{comp_id}"> F[{compress_ids(ids)}] </C>\n')
            comp_id += 1
        f.write('            <!-- fluid -->\n')
        if mesh.hexes:
            f.write(f'            <C ID="{fluid_comp_id}"> H[0-{len(mesh.hexes)-1}] </C>\n')
        else:
            f.write(f'            <C ID="{fluid_comp_id}"> H[] </C>\n')
        f.write('        </COMPOSITE>\n')

        f.write(f'        <DOMAIN> C[{fluid_comp_id}] </DOMAIN>\n')
        f.write('    </GEOMETRY>\n')

        if write_expansions:
            f.write('    <EXPANSIONS>\n')
            f.write(
                f'        <E COMPOSITE="C[{fluid_comp_id}]" '
                f'NUMMODES="{num_modes},{num_modes},{num_modes}" '
                f'BASISTYPE="Modified_A,Modified_A,Modified_A" '
                f'POINTSTYPE="GaussLobattoLegendre,GaussLobattoLegendre,GaussLobattoLegendre" '
                f'NUMPOINTS="{num_points},{num_points},{num_points}" '
                f'FIELDS="{fields}" />\n'
            )
            f.write('    </EXPANSIONS>\n')

        f.write('</NEKTAR>\n')

    print(f"Wrote {path}")
    print(f"  vertices: {len(mesh.vertices)}")
    print(f"  edges:    {len(mesh.edges)}")
    print(f"  faces:    {len(mesh.faces)}")
    print(f"  hexes:    {len(mesh.hexes)}")
    print(f"  external boundary faces: {sum(len(v) for v in bgroups.values())}")
    for name, ids in bgroups.items():
        print(f"    {name:8s}: {len(ids)}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert structured Tecplot DAT to Nektar++ 3D XML.")
    ap.add_argument("dat_files", nargs="+", help="Input Tecplot structured DAT files")
    ap.add_argument("-o", "--output", default="mesh3d.xml", help="Output Nektar++ XML file")
    ap.add_argument("--tol", type=float, default=1.0e-10, help="Coordinate merge tolerance")
    ap.add_argument("--open-j", action="store_true", help="Do not auto-detect/use j-periodic seam")
    ap.add_argument("--kleft", type=int, default=None, help="1-based kleft of the main block wall region")
    ap.add_argument("--kright", type=int, default=None, help="1-based kright of the main block wall region")
    ap.add_argument("--num-modes", type=int, default=4, help="NUMMODES in each direction")
    ap.add_argument("--num-points", type=int, default=5, help="NUMPOINTS in each direction")
    ap.add_argument("--fields", default="u,v,w,p", help="Nektar++ fields list")
    ap.add_argument("--no-expansions", action="store_true", help="Do not write EXPANSIONS section")
    args = ap.parse_args()

    mesh = Mesh3D(tol=args.tol)

    for dat in args.dat_files:
        for z in read_tecplot_dat(dat):
            build_zone(mesh, z, closed_j_auto=not args.open_j, force_open_j=args.open_j, kleft=args.kleft, kright=args.kright)

    write_xml(
        args.output,
        mesh,
        num_modes=args.num_modes,
        num_points=args.num_points,
        fields=args.fields,
        write_expansions=not args.no_expansions,
    )


if __name__ == "__main__":
    main()
