"""bam-ci: pre-processamento do BAM 1.2 subsazonal em Python.

Substitui a cadeia GrADS/Fortran que hoje gera as analises de entrada do modelo.
Estado atual (0.1.0): fundacao — leitura de descritores/binarios GrADS e
interpolacao vertical pressao -> niveis hibridos. O restante do pipeline entra
quando o codigo original estiver disponivel (jaci:/p/projetos/ioper/models/
BAM-1.2-subsaz).
"""

from .ctl import Ctl, CtlError, parse_ctl, read_ctl
from .grads import GradsFile
from .spectral import analysis as spectral_analysis
from .spectral import synthesis as spectral_synthesis
from .vertical import interp_logp, new_ps, new_sigma, vert_sigma_inter

__version__ = "0.1.0"

__all__ = [
    "Ctl", "CtlError", "parse_ctl", "read_ctl", "GradsFile",
    "interp_logp", "new_ps", "new_sigma", "vert_sigma_inter",
    "spectral_analysis", "spectral_synthesis", "__version__",
]
