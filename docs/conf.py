"""Sphinx configuration for ccc-nulls."""
import ccc_nulls

project = "ccc-nulls"
author = "ccc-nulls contributors"
copyright = "2026, ccc-nulls contributors"
release = version = ccc_nulls.__version__

extensions = ["myst_nb", "sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.viewcode",
              "sphinx.ext.intersphinx", "sphinx_copybutton"]
myst_enable_extensions = ["dollarmath", "colon_fence"]
nb_execution_mode = "off"                      # notebooks are committed executed
autodoc_typehints = "description"
autodoc_member_order = "bysource"
napoleon_numpy_docstring = True
intersphinx_mapping = {"python": ("https://docs.python.org/3", None),
                       "anndata": ("https://anndata.readthedocs.io/en/stable/", None)}
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
suppress_warnings = ["myst.header", "intersphinx"]

html_theme = "furo"
html_title = f"ccc-nulls {version}"
html_static_path = ["_static"]
