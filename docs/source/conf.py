# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

project   = 'HexMaps'
copyright = '2026, Jakob den Brok & Lukas Neumann'
author    = 'Jakob den Brok & Lukas Neumann'

# Take the version from pyproject.toml so the docs never drift from the package.
import re as _re
from pathlib import Path as _Path

_pyproject = (_Path(__file__).resolve().parents[2] / 'pyproject.toml').read_text()
release = _re.search(r'(?m)^version\s*=\s*"([^"]+)"', _pyproject).group(1)
version = release

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
