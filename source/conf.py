# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import shutil
from pathlib import Path

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'LASER: Long-read Allele-Specific Expression Resolver'
copyright = '2026, Hannah Macpherson'
author = 'Hannah Macpherson'
release = '1.1'

# -- Tutorial notebook -------------------------------------------------------
# The tutorial notebook lives in tutorial/scripts/, outside the Sphinx source
# folder, so copy it in at build time. Its figures are rendered from the
# outputs saved in the notebook, so no PNGs need committing to git.
# The copy (source/tutorial_notebook.ipynb) is gitignored.

_here = Path(__file__).resolve().parent
_nb_src = _here.parent / 'tutorial' / 'scripts' / '0_tutorial.ipynb'
_nb_dst = _here / 'tutorial_notebook.ipynb'
if _nb_src.exists():
    shutil.copyfile(_nb_src, _nb_dst)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['nbsphinx', 'sphinxcontrib.mermaid']

# Height of rendered Mermaid diagrams. The extension's 500px default squashes
# the tall pipeline flowchart on the index page until the text is unreadable.
mermaid_height = '1000px'

# Never run the notebook during the docs build: it needs the finished
# tutorial run, so show the outputs already saved in the notebook instead.
nbsphinx_execute = 'never'

templates_path = ['_templates']
exclude_patterns = ['**.ipynb_checkpoints']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    # Sidebar header (title, version and search box): pink to match the logo
    'style_nav_header_background': '#CD547D',
}
html_static_path = ['_static']
# Colour scheme (plum / pink / dusty rose / blush / charcoal)
html_css_files = ['custom.css']
