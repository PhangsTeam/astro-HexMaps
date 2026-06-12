"""
PyStructure Pipeline: homogenize and analyze multi-wavelength astronomical datasets.
"""

__author__ = "J. den Brok & L. Neumann"
__version__ = "4.3.0"
__email__ = "jadenbrok@mpia.de & lukas.neumann@eso.org"
__credits__ = ["M. Jimenez-Donaire", "E. Rosolowsky", "A. Leroy", "I. Beslic"]

from pystructurePipeline.handlerPipeline import PipelineHandler
from pystructurePipeline.init_workdir import init_workdir

__all__ = ["PipelineHandler", "init_workdir"]
