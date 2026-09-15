#!/usr/bin/env python3
"""Build folded fsaverage6 pial surfaces with projected AAL3 identifiers.

This script generates *visual geometry only*. It never reads or changes gene
expression, GWAS, enrichment, null-model, or validation statistics.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.datasets import fetch_surf_fsaverage
from scipy.ndimage import gaussian_filter, zoom
from scipy.interpolate import interpn
from skimage.measure import marching_cubes


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "data/atlases/aal_3v2/AAL3/AAL3v1.nii.gz"
REGIONS = ROOT / "data/processed/region_metadata.csv"
OUTPUT = ROOT / "web/public/data/anatomical_brain.json"
RETAINED = ROOT / "data/web/anatomical_brain.json"
REGFUSION_ASSETS = ROOT / "data/atlases/regfusion_fsaverage6"
REGFUSION_CACHE = ROOT / "data/intermediate/neuromaps"


def web_coordinates(ijk: np.ndarray, affine: np.ndarray) -> np.ndarray:
    mni = nib.affines.apply_affine(affine, ijk)
    return np.column_stack((mni[:, 0], mni[:, 2], mni[:, 1])).astype(np.float32)


def rounded_mesh(vertices: np.ndarray, faces: np.ndarray) -> dict[str, list[float] | list[int]]:
    if len(vertices) < 3 or len(faces) < 1:
        raise ValueError("An anatomical mesh has no triangles")
    return {
        "positions": np.round(vertices.astype(np.float32).reshape(-1), 2).tolist(),
        "indices": faces.astype(np.int32).reshape(-1).tolist(),
    }


def compact_faces(vertices: np.ndarray, faces: np.ndarray) -> dict[str, list[float] | list[int]]:
    used, inverse = np.unique(faces.reshape(-1), return_inverse=True)
    return rounded_mesh(vertices[used], inverse.reshape(-1, 3))


def atlas_region_surface(atlas: np.ndarray, affine: np.ndarray, region_id: int) -> dict[str, list[float] | list[int]]:
    occupied = np.argwhere(atlas == region_id)
    if len(occupied) < 2:
        raise ValueError(f"Atlas label {region_id} has fewer than two voxels")
    lo = np.maximum(occupied.min(axis=0) - 2, 0)
    hi = np.minimum(occupied.max(axis=0) + 3, atlas.shape)
    source = (atlas[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] == region_id).astype(np.float32)
    # Interpolate the label occupancy at 1 mm solely for a smooth visual shell.
    # The original integer atlas remains the authoritative identity mask.
    refined = zoom(source, 2, order=1)
    refined = gaussian_filter(refined, sigma=0.6)
    vertices, faces, _, _ = marching_cubes(refined, level=0.5, step_size=1, allow_degenerate=False)
    ijk = lo + vertices / 2.0
    return rounded_mesh(web_coordinates(ijk, affine), faces)


def registration_fusion_labels(atlas_image: nib.Nifti1Image) -> tuple[np.ndarray, np.ndarray]:
    """Project AAL3 integer labels to fsaverage5 using official RF coordinates.

    The bundled 41k coordinates were acquired using neuromaps 0.0.7. On a clean
    checkout without them, its documented `mni152_to_fsaverage(...,
    method='nearest')` API downloads the official registration-fusion archive;
    the two 41k coordinate files are then copied to stable project assets.
    """
    REGFUSION_ASSETS.mkdir(parents=True, exist_ok=True)
    coordinate_files = [REGFUSION_ASSETS / f"tpl-MNI152_space-fsaverage_den-41k_hemi-{hemi}_regfusion.txt" for hemi in ("L", "R")]
    if not all(path.is_file() for path in coordinate_files):
        os.environ["NEUROMAPS_DATA"] = str(REGFUSION_CACHE)
        from neuromaps.transforms import mni152_to_fsaverage

        mni152_to_fsaverage(atlas_image, "41k", method="nearest")
        for hemi, target in zip(("L", "R"), coordinate_files):
            source = REGFUSION_CACHE / "atlases/regfusion" / target.name
            if not source.is_file():
                raise FileNotFoundError(f"Official registration-fusion file missing: {source}")
            shutil.copy2(source, target)
    atlas = np.asarray(atlas_image.dataobj, dtype=np.int16)
    grid = tuple(range(dim) for dim in atlas.shape)
    inverse = np.linalg.inv(atlas_image.affine)
    outputs = []
    for path in coordinate_files:
        ras = np.loadtxt(path)
        ijk = nib.affines.apply_affine(inverse, ras)
        projected = interpn(grid, atlas, ijk, method="nearest", bounds_error=False, fill_value=0)
        outputs.append(np.rint(projected).astype(np.int16))
    return outputs[0], outputs[1]


def build() -> dict[str, object]:
    atlas_image = nib.load(ATLAS)
    atlas = np.asarray(atlas_image.dataobj, dtype=np.int16)
    metadata = pd.read_csv(REGIONS)
    retained = metadata.loc[metadata.retained_in_main_matrix.astype(bool)]
    retained_ids = set(retained.region_id.astype(int))
    cortical_ids = set(retained.loc[retained.broad_system.eq("cortex"), "region_id"].astype(int))

    # fsaverage6 provides actual reconstructed pial gyri and sulci. Its
    # vertex-to-MNI label correspondence comes from registration fusion, not a
    # naive nearest MNI polygon or an unrelated generic brain asset.
    pial = fetch_surf_fsaverage(mesh="fsaverage6", data_dir=ROOT / "data/intermediate/nilearn")
    projected_left, projected_right = registration_fusion_labels(atlas_image)
    surface_vertices = []
    surface_faces = []
    vertex_labels = []
    for hemi, projected in (("L", projected_left), ("R", projected_right)):
        gifti = nib.load(getattr(pial, f"pial_{'left' if hemi == 'L' else 'right'}"))
        points = np.asarray(gifti.darrays[0].data, dtype=np.float32)
        triangles = np.asarray(gifti.darrays[1].data, dtype=np.int32)
        if len(projected) != len(points):
            raise ValueError(f"Registration-fusion labels do not match fsaverage6 {hemi} vertices")
        offset = sum(len(part) for part in surface_vertices)
        surface_vertices.append(points)
        surface_faces.append(triangles + offset)
        vertex_labels.append(projected)
    vertices = np.vstack(surface_vertices)
    faces = np.vstack(surface_faces)
    vertex_labels = np.concatenate(vertex_labels)
    web = np.column_stack((vertices[:, 0], vertices[:, 2], vertices[:, 1])).astype(np.float32)

    # A triangle is assigned the majority of its three pre-projected AAL3
    # *vertex* labels. Three-way conflicts and background remain uncolored.
    triplets = vertex_labels[faces]
    labels = np.where(triplets[:, 0] == triplets[:, 1], triplets[:, 0],
                      np.where(triplets[:, 0] == triplets[:, 2], triplets[:, 0],
                               np.where(triplets[:, 1] == triplets[:, 2], triplets[:, 1], 0)))
    labels = np.where(np.isin(labels, list(cortical_ids)), labels, 0).astype(np.int16)

    meshes = []
    mapped_faces = 0
    fallback_ids = []
    for row in retained.itertuples():
        region_id = int(row.region_id)
        face_index = np.where(labels == region_id)[0] if region_id in cortical_ids else np.empty(0, dtype=np.int64)
        if len(face_index) >= 12:
            mesh = compact_faces(web, faces[face_index])
            method = "fsaverage6 pial surface; nearest AAL3 registration-fusion label"
            mapped_faces += len(face_index)
        else:
            mesh = atlas_region_surface(atlas, atlas_image.affine, region_id)
            method = "AAL3 label isosurface; 1 mm occupancy interpolation"
            fallback_ids.append(region_id)
        meshes.append({
            "region_id": region_id, "source": method,
            "voxel_count": int(np.count_nonzero(atlas == region_id)), **mesh,
        })

    return {
        "schema_version": "2.0.0",
        "atlas": "AAL3v1",
        "coordinate_system": "WebGL axes (RAS x, RAS z, RAS y), millimetres; cortex uses fsaverage6 surface coordinates, atlas isosurfaces use MNI152 coordinates",
        "geometry_source": "Nilearn fsaverage6 pial cortical surface; neuromaps 0.0.7 MNI152 registration fusion; project AAL3v1 2 mm labels",
        "label_rule": "nearest MNI-to-fsaverage projection of original AAL3v1 integer labels; triangle majority vote across its projected vertex IDs; original atlas integer ID is statistical join key",
        "region_count": len(meshes),
        "anatomy": rounded_mesh(web, faces),
        "regions": meshes,
        "quality": {
            "anatomy_vertices": len(vertices), "anatomy_faces": len(faces),
            "cortical_faces_with_retained_aal3_labels": mapped_faces,
            "atlas_isosurface_region_ids": fallback_ids,
            "unassigned_anatomy_faces": int(np.count_nonzero(~np.isin(labels, list(retained_ids)))),
        },
    }


def main() -> None:
    geometry = build()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    RETAINED.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(geometry, separators=(",", ":")) + "\n"
    OUTPUT.write_text(encoded, encoding="utf-8")
    RETAINED.write_text(encoded, encoding="utf-8")
    print(f"Built {geometry['quality']['anatomy_faces']:,} anatomical faces and {geometry['region_count']} AAL3 mappings")
    print(f"Web geometry: {OUTPUT.stat().st_size / 1024 / 1024:.2f} MiB")


if __name__ == "__main__":
    main()
