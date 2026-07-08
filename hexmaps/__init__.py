"""
HexMaps Pipeline: homogenize and analyze multi-wavelength astronomical datasets.
"""

__author__ = "J. den Brok & L. Neumann"
__version__ = "5.0.0"
__email__ = "jadenbrok@mpia.de & lukas.neumann@eso.org"
__credits__ = ["M. Jimenez-Donaire", "E. Rosolowsky", "A. Leroy", "I. Beslic"]

from hexmaps.handler_pipeline import PipelineHandler
from hexmaps.init_workdir import init_workdir
from hexmaps.hexmaps_analysis import HexMapsAnalysis

__all__ = ["PipelineHandler", "init_workdir", "HexMapsAnalysis"]
