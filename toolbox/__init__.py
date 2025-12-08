# SPDX-License-Identifier: LGPL-3.0-or-later
r"""Toolbox for computational chemistry and materials science.

This package provides tools for:
- Calculator interfaces (CP2K, LAMMPS, Deep Potential)
- Input/output operations for various file formats
- Plotting and visualization utilities
- Mathematical and optimization tools
- Machine learning utilities

Examples
--------
>>> import toolbox
>>> toolbox.test()
 _       __     __                             __           __  __                 __      __
| |     / /__  / /________  ____ ___  ___     / /_____     / /_/ /_  ___     _____/ /_  __/ /_
| | /| / / _ \/ / ___/ __ \/ __ `__ \/ _ \   / __/ __ \   / __/ __ \/ _ \   / ___/ / / / / __ \\
| |/ |/ /  __/ / /__/ /_/ / / / / / /  __/  / /_/ /_/ /  / /_/ / / /  __/  / /__/ / /_/ / /_/ /
|__/|__/\___/_/\___/\____/_/ /_/ /_/\___/   \__/\____/   \__/_/ /_/\___/   \___/_/\__,_/_.___/
"""
import json
import os

fname = os.path.join(__path__[0], "../config.json")
if os.path.exists(fname):
    with open(fname) as f:
        CONFIGS = json.loads(f.read())
else:
    CONFIGS = {}


def test():
    """Print toolbox ASCII art logo."""
    print(
        " _       __     __                             __           __  __                 __      __  "
    )
    print(
        "| |     / /__  / /________  ____ ___  ___     / /_____     / /_/ /_  ___     _____/ /_  __/ /_ "
    )
    print(
        "| | /| / / _ \/ / ___/ __ \/ __ `__ \/ _ \   / __/ __ \   / __/ __ \/ _ \   / ___/ / / / / __ \\"
    )
    print(
        "| |/ |/ /  __/ / /__/ /_/ / / / / / /  __/  / /_/ /_/ /  / /_/ / / /  __/  / /__/ / /_/ / /_/ /"
    )
    print(
        "|__/|__/\___/_/\___/\____/_/ /_/ /_/\___/   \__/\____/   \__/_/ /_/\___/   \___/_/\__,_/_.___/ "
    )
