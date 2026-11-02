"""Self-contained image results with explicit provenance and table storage."""

import json
import h5py
import numpy as np
import pandas as pd


def write_json(h5, name, value):
    h5.create_dataset(
        name,
        data=json.dumps(value, default=str, allow_nan=False),
        dtype=h5py.string_dtype("utf-8"),
    )


def write_table(h5, name, frame):
    group = h5.create_group(name)
    group.attrs["columns_json"] = json.dumps(list(frame.columns))
    for column in frame:
        series = frame[column]
        if pd.api.types.is_numeric_dtype(series.dtype) or pd.api.types.is_bool_dtype(
            series.dtype
        ):
            group.create_dataset(column, data=series.to_numpy())
        else:
            group.create_dataset(
                column,
                data=series.fillna("").astype(str).to_numpy(),
                dtype=h5py.string_dtype("utf-8"),
            )


def read_table(h5, name):
    group = h5[name]
    columns = json.loads(group.attrs["columns_json"])
    values = {}
    for column in columns:
        dataset = group[column]
        values[column] = (
            dataset.asstr()[:] if h5py.check_string_dtype(dataset.dtype) else dataset[:]
        )
    frame = pd.DataFrame(values)
    # Older tables can store numeric identifiers as UTF-8 strings when their
    # pandas columns originally had object dtype. Normalize identifiers only;
    # marker names, status strings and other text retain their original type.
    for column in (
        "soma_id",
        "primary_nucleus_id",
        "proposed_nucleus_id",
        "punctum_id",
        "owner_soma_id",
    ):
        if column not in frame or pd.api.types.is_numeric_dtype(frame[column]):
            continue
        numeric = pd.to_numeric(frame[column].replace("", np.nan), errors="raise")
        present = numeric.dropna().to_numpy()
        if not np.isfinite(present).all() or np.any(present != np.floor(present)):
            raise ValueError("Non-integer identifier in " + name + "/" + column)
        frame[column] = numeric
    return frame
