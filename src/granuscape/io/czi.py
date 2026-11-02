"""Native CZI planes and unmodified acquisition metadata."""

import hashlib
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import czifile


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_metadata(xml, channel_mapping):
    root = ET.fromstring(xml)
    # Namespaces are removed only from the parsing copy. Original XML is retained.
    for element in root.iter():
        element.tag = element.tag.rsplit("}", 1)[-1]
    channels = root.findall(".//Information/Image/Dimensions/Channels/Channel")
    names = [c.get("Name", "") for c in channels]
    indices = {}
    for semantic, name in channel_mapping.items():
        if names.count(name) != 1:
            raise ValueError(
                f"CZI channel must occur exactly once: {name}; found {names}"
            )
        indices[semantic] = names.index(name)
    distances = {
        d.get("Id"): float(d.findtext("Value")) * 1e6
        for d in root.findall(".//Scaling/Items/Distance")
        if d.findtext("Value")
    }
    if not {"X", "Y"} <= distances.keys():
        raise ValueError("Physical XY calibration is missing")
    if distances["X"] <= 0 or not np.isclose(distances["X"], distances["Y"], rtol=1e-6):
        raise ValueError("Only positive isotropic XY calibration is supported")
    return indices, distances["X"]


def read_czi(path, channel_mapping):
    with czifile.CziFile(path) as czi:
        if hasattr(czi, "scenes"):
            if len(czi.scenes) != 1:
                raise ValueError("Only single-scene CZI is supported")
            image = czi.scenes[0]
        else:
            image = czi
        data, axes, xml = image.asarray(), image.axes, czi.metadata()
        if isinstance(xml, bytes):
            xml = xml.decode("utf-8")
        for axis in list(axes):
            if axis not in "CYX":
                index = axes.index(axis)
                if data.shape[index] != 1:
                    raise ValueError(
                        f"Non-singleton {axis} axis: no implicit projection is performed"
                    )
                data = np.take(data, 0, axis=index)
                axes = axes[:index] + axes[index + 1 :]
        if set(axes) != set("CYX") or len(axes) != 3:
            raise ValueError(f"Unsupported CZI axes: {axes}")
        data = np.transpose(data, [axes.index(axis) for axis in "CYX"])
        if data.dtype != np.uint16:
            raise ValueError("Native uint16 fluorescence planes are required")
        indices, pixel_size = parse_metadata(xml, channel_mapping)
        if max(indices.values()) >= len(data):
            raise ValueError("Channel metadata and array dimensions differ")
    return {name: data[i] for name, i in indices.items()}, pixel_size, xml
