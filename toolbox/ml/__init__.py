# SPDX-License-Identifier: LGPL-3.0-or-later
"""Machine learning utilities for computational chemistry.

This module provides tools for:
- TensorFlow graph visualization and analysis
- Deep Potential model integration

Examples
--------
>>> from toolbox.ml import Graph
>>> graph = Graph("model.pb")
>>> graph.run()  # Launch TensorBoard visualization
"""
from .tf_graph import Graph

__all__ = ["Graph"]
