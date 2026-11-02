"""Best-effort hardware inventory; unavailable fields remain explicit."""

import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
from pathlib import Path


def query(command, timeout=30):
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False
        )
        return {
            "returncode": result.returncode,
            "stdout": result.stdout.strip(),
            "stderr": result.stderr.strip(),
        }
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"error": str(error)}


def inventory(config, output):
    packages = {}
    for name in (
        "numpy",
        "scipy",
        "pandas",
        "h5py",
        "scikit-image",
        "scikit-learn",
        "czifile",
        "Pillow",
    ):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = "unavailable"
    info = {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cpu": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "packages": packages,
        "gpu_requested": config.get("gpu", True),
        "power_protocol": "Use AC power; prevent sleep; close unrelated compute tasks. Power/thermal conditions are not automatically controlled.",
        "ram_measurement": "Installed capacity only, not peak process memory",
    }
    if os.name == "nt":
        command = "$c=Get-CimInstance Win32_Processor; $s=Get-CimInstance Win32_ComputerSystem; [pscustomobject]@{CPU=($c.Name -join '; '); PhysicalCores=($c.NumberOfCores | Measure-Object -Sum).Sum; LogicalCores=($c.NumberOfLogicalProcessors | Measure-Object -Sum).Sum; RAMBytes=$s.TotalPhysicalMemory} | ConvertTo-Json -Compress"
        info["windows_hardware"] = query(
            ["powershell", "-NoProfile", "-Command", command]
        )
    else:
        try:
            info["installed_ram_bytes"] = os.sysconf("SC_PAGE_SIZE") * os.sysconf(
                "SC_PHYS_PAGES"
            )
        except (ValueError, OSError, AttributeError):
            info["installed_ram_bytes"] = None
    nvidia = shutil.which("nvidia-smi")
    info["nvidia"] = (
        query(
            [
                nvidia,
                "--query-gpu=name,driver_version,memory.total",
                "--format=csv,noheader",
            ]
        )
        if nvidia
        else {"error": "nvidia-smi unavailable; use recorded Cellpose device metadata"}
    )
    info["output_disk"] = {
        "path": str(output),
        "free_bytes": shutil.disk_usage(output).free,
    }
    info["thread_environment"] = {
        k: os.environ.get(k)
        for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")
    }
    return info
