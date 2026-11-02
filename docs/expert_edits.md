# Image-bound expert edits

Set `edits_json` in the YAML configuration to an expert-edits JSON file. Each entry is bound to the SHA-256 of the original CZI; a mismatch stops that image. Image-specific decisions are not embedded in numerical modules.

```json
{
  "images": [
    {
      "image_id": "example_image",
      "source_sha256": "replace_with_full_sha256_of_original_czi",
      "merges": [{"source_id": 12, "target_id": 11}],
      "mcl_boundary_xy_px": [[10, 0], [12, 100], [14, 200]]
    }
  ]
}
```

`merges` and `mcl_boundary_xy_px` are independently optional. Coordinates are native image pixels in x,y order. The illustrated coordinates are a format example, not an anatomically approved boundary. Curves must be open and free of self-intersections. Merges must be connected, with both labels present and neither touching the frame. Original, merged and refined label arrays remain separate in the result HDF5. The exact expert edit record is included in provenance.
