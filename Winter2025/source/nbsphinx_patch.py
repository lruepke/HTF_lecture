# _extensions/nbsphinx_patch.py
from nbconvert.exporters.html import HTMLExporter   # nbsphinx uses HTMLExporter
import sys

def _hide_cells_inplace(nb):
    if not getattr(nb, "cells", None):
        return
    print("nbsphinx_patch: scanning for 'hide'…")
    for i, cell in enumerate(nb.cells):
        src = (cell.get("source") or "")
        if not src:
            continue
        first = src.split("\n", 1)[0]
        if first and "hide" in first:
            print(f"nbsphinx_patch: hiding cell #{i}")
            msg = f"{first}: the code in this cell is hidden by the author"
            cell["source"] = msg
            # optionally also clear outputs:
            # cell["outputs"] = []
            # cell["execution_count"] = None

# Patch *the* method HTMLExporter calls:
_ORIG = HTMLExporter.from_notebook_node

def _patched_from_notebook_node(self, nb, resources=None, **kw):
    _hide_cells_inplace(nb)
    return _ORIG(self, nb, resources=resources, **kw)

HTMLExporter.from_notebook_node = _patched_from_notebook_node
print("nbsphinx_patch: HTMLExporter.from_notebook_node patched", file=sys.stderr)
