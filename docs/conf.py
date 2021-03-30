# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import sys
from os.path import dirname, join

SRC_PATH = join(dirname(dirname(__file__)), "geckopy")
sys.path.insert(0, SRC_PATH)


# -- Project information -----------------------------------------------------

project = 'geckopy'
# copyright = '
# author = ''

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'autoapi.extension',
]

autoapi_dirs = [SRC_PATH]

intersphinx_mapping = {'cobrapy': ('https://cobrapy.readthedocs.io/en/latest/', None)}
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'
html_favicon = '_static/favicon.ico'
html_static_path = ['_static']

html_theme_options = {
    # 'logo': 'geckopy_lojo.png',
    'logo': 'logo.png',
    'touch_icon': 'logo.png',
    'github_banner': True,
    'github_user': 'ginkgobioworks',
    'github_repo': 'geckopy',
    'page_width': '1200px'
}
mathjax_path = (
    "https://cdn.mathjax.org/mathjax/latest/"
    "MathJax.js?config=TeX-AMS-MML_HTMLorMML"
)

# -- Options for LaTeX output --------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    "papersize": "a4paper",
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    "preamble": r"\usepackage{amsmath,amssymb}",
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
