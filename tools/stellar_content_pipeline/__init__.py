"""Stellar content pipeline.

This package is meant to live inside the repo under tools/ so it can be used
during development without external dependencies.

Run:
  python -m tools.stellar_content_pipeline --input tools/examples/universe_state.json --output out/content
"""

from .cli import main  # re-export convenience
