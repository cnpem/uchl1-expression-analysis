"""Process-wide preload of the workbench's heavy imports.

`marimo run` serves every browser tab as a session thread in one process. When two sessions start
together they can import pandas, IPython or anywidget concurrently, and marimo's import hooks then
see half-initialised modules. `preload()` imports them once, under one lock, before any notebook
code uses them.
"""

from __future__ import annotations

import threading

IMPORT_LOCK = threading.Lock()


def preload() -> None:
    with IMPORT_LOCK:
        try:
            import IPython.display  # noqa: F401
        except ImportError:
            pass
        import pandas  # noqa: F401

        try:
            import anywidget  # noqa: F401
        except ImportError:
            pass
        import altair  # noqa: F401
