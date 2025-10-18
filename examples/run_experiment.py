# -*- coding: utf-8 -*-
"""
Example: Run McGhee–von Hippel analysis across multiple titrations (replicates).
"""

from pathlib import Path
from mvh_analysis.Experiment import Experiment

config_paths = [
    Path("examples/config_example_rep1.toml"),
    Path("examples/config_example_rep2.toml"),
]

exp = Experiment(config_paths)
exp.initialize_all()
exp.run_all()
summary = exp.summarize()
averages = exp.average_results()
exp.plot_all()

print("\nSummary table:")
print(summary)
