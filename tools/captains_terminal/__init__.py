"""Captain's Terminal: an interactive, text-mode gameplay loop.

This is a *developer tool* that sits next to the offline content pipeline.
It exists to make the project playable *today* even when the in-engine UI/
mission systems are still under construction.

Features:
- Procedural GalNet feed + mission board (driven by stellar_content_pipeline)
- Accept/complete/fail missions
- Travel between systems with risk-driven encounters
- Lightweight progression (credits, hull, faction reputation)

Run:
    python -m tools.captains_terminal --help
"""

from __future__ import annotations

__all__ = ["__version__"]

__version__ = "0.1.0"
