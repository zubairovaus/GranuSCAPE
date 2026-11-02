"""Native-grid QC views with display-only normalization and readable cell IDs."""

from functools import lru_cache
from io import BytesIO
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi
from skimage import exposure
from skimage.segmentation import find_boundaries

from .settings import resolve_options

COLORS = {
    "soma": "#FFD65A",
    "nuclei": "#80CBFF",
    "g3bp": "#FF655F",
    "eif3": "#46F19A",
    "intersection": "#FFFFFF",
}
VIEW_NAMES = {
    "A": "A_overlay.png",
    "B": "B_masks.png",
    "C": "C_cells.png",
    "D": "D_puncta.png",
}
TITLES = {
    "A": "RGB + bright contours",
    "B": "Black + translucent masks",
    "C": "RGB + soma / nuclei",
    "D": "Black + puncta only",
}


@lru_cache(maxsize=16)
def font(size):
    for name in ("DejaVuSans.ttf", "Arial.ttf", "C:/Windows/Fonts/arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def rgb_color(value):
    return (
        np.array([int(value[i : i + 2], 16) for i in (1, 3, 5)], dtype=np.float32) / 255
    )


def blend(image, mask, color, alpha):
    image[mask] = (1 - alpha) * image[mask] + alpha * rgb_color(color)


def composite(raw, options=None):
    """Normalize full channels for display, never for masks or measurements."""
    options = resolve_options(options)
    channels = []
    for name in ("rhoda", "af488", "dapi"):
        array = np.asarray(raw[name], dtype=np.float32)
        if array.ndim != 2 or not np.isfinite(array).all():
            raise ValueError("QC requires finite two-dimensional raw channels")
        low, high = np.percentile(array, [1, 99.7])
        if high <= low:
            channels.append(np.zeros_like(array))
            continue
        normalized = np.clip((array - low) / max(high - low, 1), 0, 1)
        if options["normalization"] == "clahe":
            normalized = exposure.equalize_adapthist(
                normalized,
                kernel_size=tuple(options["clahe_kernel_px"]),
                clip_limit=options["clahe_clip_limit"],
                nbins=256,
            ).astype(np.float32)
            normalized **= options["clahe_gamma"]
        channels.append(normalized)
    if len({c.shape for c in channels}) != 1:
        raise ValueError("QC raw channel dimensions differ")
    return np.stack(channels, axis=-1)


def render_views(raw, soma, nuclei, granules, roi, options):
    shape = soma.shape
    for array in (nuclei, granules["g3bp"], granules["eif3"], roi):
        if array.shape != shape:
            raise ValueError("QC mask dimensions differ")
    wanted = set(options["views"])
    if options["include_cell_ids"]:
        wanted.add("A")
    if not wanted:
        return
    rgb = composite(raw, options) if wanted & {"A", "C"} else None
    edges = {
        name: find_boundaries(array, mode="inner")
        for name, array in (
            ("soma", soma),
            ("nuclei", nuclei),
            ("g3bp", granules["g3bp"]),
            ("eif3", granules["eif3"]),
        )
    }
    yy, xx = np.indices(shape)
    border = find_boundaries(roi, mode="inner")
    dilation = (options["mcl_line_width_px"] - 1) // 2
    if dilation:
        border = ndi.binary_dilation(border, iterations=dilation)
    border &= (xx + yy) % 16 < 10
    overlap = (granules["g3bp"] > 0) & (granules["eif3"] > 0)
    for view in "ABCD":
        if view not in wanted:
            continue
        image = (
            rgb * options["background_brightness"]
            if view in ("A", "C")
            else np.zeros(shape + (3,), dtype=np.float32)
        )
        included = (
            ["soma", "nuclei"]
            if view == "C"
            else ["g3bp", "eif3"] if view == "D" else ["soma", "nuclei", "g3bp", "eif3"]
        )
        if view == "B":
            blend(image, soma > 0, COLORS["soma"], 0.16)
            blend(image, nuclei > 0, COLORS["nuclei"], 0.48)
        if view in ("B", "D"):
            for marker in ("g3bp", "eif3"):
                blend(image, granules[marker] > 0, COLORS[marker], 0.86)
        if view in ("A", "C"):
            keyline = np.logical_or.reduce(
                [ndi.binary_dilation(edges[key], iterations=1) for key in included]
            )
            blend(image, keyline, "#000000", 0.82)
        for key in included:
            edge = (
                edges[key]
                if key == "nuclei"
                else ndi.binary_dilation(edges[key], iterations=1)
            )
            blend(image, edge, COLORS[key], 0.97)
        if view != "C":
            blend(image, overlap, COLORS["intersection"], 1.0)
        blend(image, border, options["mcl_color"], 1.0)
        yield view, Image.fromarray(
            np.rint(np.clip(image, 0, 1) * 255).astype(np.uint8)
        )


def decorate(image, view, curve, cells, pixel_size, options, numbered=False):
    image = image.copy()
    draw = ImageDraw.Draw(image)
    width, height = image.size
    stroke = options["mcl_line_width_px"]
    points = [(float(x), float(y)) for x, y in curve]
    if len(points) >= 2:
        draw.line(points, fill="#101820", width=stroke + 2)
        draw.line(points, fill=options["mcl_color"], width=stroke)
    if numbered:
        required = {
            "soma_id",
            "soma_centroid_x_px",
            "soma_centroid_y_px",
            "analysis_ROI_fraction",
        }
        if not required <= set(cells.columns):
            raise ValueError("Cell table lacks columns required for QC identifiers")
        for row in cells.itertuples():
            if row.analysis_ROI_fraction < 0.1:
                continue
            x, y = row.soma_centroid_x_px, row.soma_centroid_y_px
            if not (
                np.isfinite(x) and np.isfinite(y) and 0 <= x < width and 0 <= y < height
            ):
                continue
            label = str(int(row.soma_id))
            bounds = draw.textbbox((x, y), label, font=font(14), anchor="mm")
            draw.rectangle(
                (bounds[0] - 2, bounds[1] - 1, bounds[2] + 2, bounds[3] + 1),
                fill="#101820",
            )
            draw.text((x, y), label, font=font(14), anchor="mm", fill="white")
    # Draw calibration last so labels cannot obscure it.
    available_um = max(1, width - 32) * pixel_size
    length_um = next(
        (v for v in (50, 20, 10, 5, 2, 1) if v < available_um * 0.45),
        available_um * 0.25,
    )
    length = round(length_um / pixel_size)
    pad = min(max(12, round(width * 0.025)), max(2, min(width, height) // 5))
    end, y = width - pad, height - pad
    draw.line((end - length, y, end, y), fill="black", width=10)
    draw.line((end - length, y, end, y), fill="white", width=5)
    draw.text(
        (end - length / 2, y - 12),
        f"{length_um:g} µm",
        font=font(max(14, round(width * 0.011))),
        anchor="ms",
        fill="white",
        stroke_width=2,
        stroke_fill="black",
    )
    keys = (
        ["soma", "nuclei"]
        if view == "C"
        else ["g3bp", "eif3", "intersection"] if view == "D" else list(COLORS)
    )
    keys += ["MCL ROI", "MCL interface"]
    names = {
        "soma": "Soma",
        "nuclei": "Nuclei",
        "g3bp": "G3BP",
        "eif3": "eIF3η",
        "intersection": "Pixel overlap",
    }
    layout = []
    x, y = 16, 52 if width > 1000 else 37
    face = font(19 if width > 1000 else 13)
    for key in keys:
        text = names.get(key, key)
        color = COLORS.get(key, options["mcl_color"])
        span = draw.textlength(text, font=face) + 36
        if x + span > width - 16:
            x, y = 16, y + 26
        layout.append((x, y, color, text))
        x += span
    header = max(112 if width > 1000 else 120, y + 28)
    canvas = Image.new("RGB", (width, height + header), "#101820")
    canvas.paste(image, (0, header))
    draw = ImageDraw.Draw(canvas)
    title = TITLES[view] + (" | cell IDs in MCL" if numbered else "")
    draw.text((16, 10), title, font=font(24 if width > 1000 else 16), fill="white")
    for x, y, color, text in layout:
        draw.rectangle((x, y + 4, x + 13, y + 15), fill=color)
        draw.text((x + 20, y), text, font=face, fill="white")
    return canvas


def save_png(image, path):
    """Publish complete PNG streams atomically, including all native image pixels."""
    stream = BytesIO()
    image.save(stream, format="PNG")
    data = stream.getvalue()
    if not data.endswith(b"\x00\x00\x00\x00IEND\xaeB`\x82"):
        raise RuntimeError("Incomplete PNG encoding")
    temporary = path.with_suffix(".png.partial")
    temporary.write_bytes(data)
    temporary.replace(path)


def save_qc(
    raw, soma, nuclei, granules, curve, roi, cells, pixel_size, folder, options=None
):
    options = resolve_options(options)
    if not np.isfinite(pixel_size) or pixel_size <= 0:
        raise ValueError("QC requires positive finite pixel calibration")
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    outputs = []
    for view, image in render_views(raw, soma, nuclei, granules, roi, options):
        if view in options["views"]:
            name = VIEW_NAMES[view]
            save_png(
                decorate(image, view, curve, cells, pixel_size, options), folder / name
            )
            outputs.append(name)
        if view == "A" and options["include_cell_ids"]:
            save_png(
                decorate(image, view, curve, cells, pixel_size, options, numbered=True),
                folder / "cell_ids.png",
            )
            outputs.append("cell_ids.png")
    settings = {
        "options": options,
        "colors": COLORS,
        "outputs": outputs,
        "native_shape_yx": list(soma.shape),
        "pixel_size_um": float(pixel_size),
        "normalization_scope": "full channels before composition; display only",
        "raw_or_masks_modified": False,
        "ids_scope": "all existing soma IDs with analysis_ROI_fraction >= 0.1; no morphotype implication",
        "MCL_style": "same color and width; solid interface, dashed ROI",
        "display_contours": "one-pixel emphasis for soma and puncta, not a change to masks or areas",
        "header_outside_image": True,
    }
    (folder / "display_settings.json").write_text(
        json.dumps(settings, indent=2), encoding="utf-8"
    )
    return outputs
