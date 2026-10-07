"""Paths that locate the reused validated modules in HEP_neurofeedback/python.

get_rr_v25.py and rr_cleaning_v25.py live there (they were written and
validated for that project) and are imported from here so the BrainBeats port
has one canonical copy of each. If they ever move, update _EXTRA_REPO_ROOTS.
"""
from __future__ import annotations

import importlib
import importlib.util
import sys
from pathlib import Path

_EXTRA_REPO_ROOTS = [
    Path(r'C:\Users\ccann\Documents\HEP_neurofeedback\python'),
]

_HERE = Path(__file__).resolve().parent


def load_module(name: str):
    """Import `name` from HEP_neurofeedback/python, falling back to normal import."""
    for root in _EXTRA_REPO_ROOTS:
        path = root / f'{name}.py'
        if path.exists():
            if name in sys.modules:
                return sys.modules[name]
            if str(root) not in sys.path:      # siblings use bare imports
                sys.path.insert(0, str(root))
            spec = importlib.util.spec_from_file_location(name, path)
            mod = importlib.util.module_from_spec(spec)
            sys.modules[name] = mod
            spec.loader.exec_module(mod)
            return mod
    return importlib.import_module(name)


def _getrr():
    m = load_module('get_rr_v25')
    return m.get_rr, m.resolve_polarity


def _cleanrr():
    m = load_module('rr_cleaning_v25')
    return m.clean_rr, m.qrs_bandpass