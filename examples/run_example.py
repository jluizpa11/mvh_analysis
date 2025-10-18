# -*- coding: utf-8 -*-
"""
Example: Run a McGhee–von Hippel fit from a configuration file.
"""

from pathlib import Path
from mvh_analysis.Titration import Titration

config_path = Path("examples/config_example.toml")

t = Titration(config_path)
t.initialize_mvh_analysis()
