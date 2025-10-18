🧬 McGhee–von Hippel (MVH) Analysis



A Python package for fitting DNA–ligand binding isotherms using the

McGhee–von Hippel model for cooperative and non-cooperative binding.



🚀 Overview



This package streamlines the quantitative analysis of ligand–DNA titrations

(absorbance, fluorescence, etc.) to determine:



Binding constant (K)



Binding site size (n)



Goodness of fit (R²)



Optional cooperativity parameter (ω) in future versions



It provides both low-level analytical functions and high-level workflow classes

(Titration, Experiment) for single and replicate analyses.



🧩 Features



🔬 Fit McGhee–von Hippel binding model to titration data



⚙️ Config-driven reproducibility via .toml configuration files



📈 Automatic generation of publication-quality plots



🧮 Easy replicate handling and averaging



🧾 CSV import/export for integration with UVTools and spectroscopy data



🧪 Quick Start

Single titration

from pathlib import Path

from mvh\_analysis.Titration import Titration



config\_path = Path("examples/config\_example.toml")

t = Titration(config\_path)

t.initialize\_mvh\_analysis()



Multiple replicates

from pathlib import Path

from mvh\_analysis.Experiment import Experiment



configs = \[Path("rep1.toml"), Path("rep2.toml"), Path("rep3.toml")]

exp = Experiment(configs)

exp.initialize\_all()

exp.run\_all()

exp.summarize()

exp.average\_results()

exp.plot\_all()



⚙️ Configuration Example (config\_example.toml)

\[paths]

filepath = "examples/example\_titration.csv"



\[experiment]

wavelength\_nm = 506

conc\_DNA\_uM = 1150.0

conc\_dye\_uM = 44.0

initial\_volume\_uL = 500.0



\[fit]

p0 = \[1e6, 2.0]

bounds\_lower = \[1e3, 0.5]

bounds\_upper = \[1e8, 10.0]

function = "mvh\_rhs"



📊 Typical Output

K = 1.26e+06 ± 1.45e+05 M⁻¹

n = 1.66 ± 0.08

R² = 0.9436



📚 Reference



McGhee, J. D.; von Hippel, P. H.

Theoretical Aspects of DNA–Protein Interactions: Cooperative and Non-Cooperative Binding of Large Ligands to a One-Dimensional Homogeneous Lattice.

J. Mol. Biol., 1974, 86, 469–489.



🧑‍💻 Author



João Luiz Petrarca de Albuquerque

University of Pennsylvania — Chenoweth Lab

📧 jluizpa11@gmail.com

