"""Shared, disk-persisted cache for biohub SAE feature lookups.

WHY: every HP3 script that needs a feature's label/category/description
(``hp3_feature_lookup.py``, ``residue_clusters_hp3.py``, and the upcoming
whole-protein scan) was hitting the network with its own in-memory
``lru_cache`` -- so the same feature id got re-fetched from scratch every time
a script reran, even across different scripts in the same session. This
module is a single on-disk cache (JSON, feature_id -> full API response) that
all of them share: a feature is fetched from the network at most once, ever,
across the whole project.

Usage:
    from hp3_feature_cache import get_feature_info
    info = get_feature_info(6113)   # dict with label/category/description/...
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import requests

CACHE_PATH = Path(__file__).resolve().parent.parent / "analysis" / "features" / "biohub_feature_cache.json"

_cache: dict[str, dict[str, Any]] | None = None


def _load() -> dict[str, dict[str, Any]]:
    global _cache
    if _cache is None:
        if CACHE_PATH.exists():
            _cache = json.loads(CACHE_PATH.read_text())
        else:
            _cache = {}
    return _cache


def _save() -> None:
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    CACHE_PATH.write_text(json.dumps(_cache, indent=0))


def get_feature_info(feature_idx: int) -> dict[str, Any]:
    """Fetch metadata for one SAE feature, disk-cached across all HP3 scripts."""
    cache = _load()
    key = str(feature_idx)
    if key in cache:
        return cache[key]
    url = f"https://biohub.ai/esm/protein/api/v1alpha1/features/{feature_idx}"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    info = response.json()
    cache[key] = info
    _save()
    return info


def cache_size() -> int:
    return len(_load())
