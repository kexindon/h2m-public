 # Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'H2M'
copyright = '2024, Kexin Dong'
author = 'Kexin Dong'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.duration",
    'sphinx.ext.autodoc',
]

templates_path = ['_templates']
exclude_patterns = [] 
import sys
import os

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static']
html_permalinks_icon = '<span>#</span>'
html_theme = 'sphinxawesome_theme'
html_theme_options = {
    "collapse_navigation": True,
    "logo_light": "figures/h2m-logo-final.png",
    "logo_dark": "figures/h2m-logo-final.png"
}

html_css_files = [
    'css/custom.css',
]