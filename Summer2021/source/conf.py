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
import sys, os
sys.path.append(os.path.abspath('_extensions'))


# -- Project information -----------------------------------------------------

project = 'Finite Element Modeling in Geodynamics'
latex_name='LectureNote'
copyright = '2021, Lars Ruepke and Zhikui Guo'
author = 'Lars Ruepke'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'nbsphinx',
    'sphinx.ext.mathjax',
    'sphinxcontrib.bibtex',
    'sphinx_inline_tabs',
    'matplotlib.sphinxext.plot_directive',
]
bibtex_bibfiles = ['refs.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '**.ipynb_checkpoints']

# config plot 
plot_html_show_source_link=True
plot_include_source=False
plot_formats=['svg','pdf']
plot_pre_code="""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
mpl.rcParams["font.family"] = 'Arial' 
mpl.rcParams["mathtext.fontset"] = 'cm'
"""

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}

html_context = {
    # Enable the "Edit in GitHub link within the header of each page.
    'display_github': True,
    # Set the following variables to generate the resulting github URL for each page.
    # Format Template: https://{{ github_host|default("github.com") }}/{{ github_user }}/{{ github_repo }}/blob/{{ github_version }}{{ conf_py_path }}{{ pagename }}{{ suffix }}
    'github_user': 'lruepke',
    'github_repo': 'HTF_lecture',
    'github_version': 'main/',
    'conf_py_path': '/source/',
    "menu_links": [
        (
            '<i class="fa fa-envelope fa-fw"></i> Prof. Lars Ruepke',
            "mailto:lruepke@geomar.de",
            
        ),
        (
            '<i class="fa fa-envelope fa-fw"></i> Dr. Zhikui Guo',
            "mailto:zguo@geomar.de",
            
        ),
        (
            '<i class="fa fa-home fa-fw"></i> Public site',
            "https://lruepke.github.io/FEM_lecture/",
        ),
        (
            '<i class="fa fa-home fa-fw"></i> University site',
            "https://lms.uni-kiel.de/url/RepositoryEntry/3664576582",
        ),
    ],
    'project':project,
    'downloads_url':'https://lruepke.github.io/FEM_lecture/downloads',
    'latex_main':  latex_name, #pdf file name
    'pdf_versions': [
        (
            'latest',
            'https://lruepke.github.io/FEM_lecture/downloads'  #pdf path
        ),
        (
            '2.0',
            '#'
        ),
        ]
}


html_logo="_static/logo/htf_picto2.png"
# some customized css style
html_css_files = ["style.css"]
# customize OpenFOAM syntax highlight
from sphinx.highlighting import lexers
from pygments_OpenFOAM.foam import OpenFOAMLexer
lexers['foam'] = OpenFOAMLexer(startinline=True)
# default language to highlight source code
highlight_language = 'foam'
pygments_style = 'xcode' # xcode,monokai,emacs,autumn,vs,solarized-dark
# numbering of cross reference
numfig = True 
math_numfig = True
math_eqref_format = 'equation ({number})'
# only valid for latex
# numfig_format = 'Figure. %s'
numfig_secnum_depth = 1
imgmath_latex = 'dvilualatex'
imgmath_image_format = 'svg'
imgmath_dvipng_args = ['-gamma', '1.5', '-bg', 'Transparent']

# -- Options for LaTeX output ---------------------------------------------
latex_engine = 'xelatex'
latex_elements = {
    'papersize': 'a4paper',
    'utf8extra': '',
    'inputenc': '',
    'babel': r'''\usepackage[english]{babel}''',
}
master_doc = 'index'
latex_documents = [
    (master_doc, latex_name+'.tex', project, author, 'manual'),
]
# ====================================================================
# new defined cite style
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.plugin import register_plugin
from collections import Counter
import re
import unicodedata

from pybtex.style.labels import BaseLabelStyle

_nonalnum_pattern = re.compile('[^A-Za-z0-9 \-]+', re.UNICODE)

def _strip_accents(s):
    return "".join(
        (c for c in unicodedata.normalize('NFD', s)
            if not unicodedata.combining(c)))

def _strip_nonalnum(parts):
    """Strip all non-alphanumerical characters from a list of strings.

    >>> print(_strip_nonalnum([u"ÅA. B. Testing 12+}[.@~_", u" 3%"]))
    AABTesting123
    """
    s = "".join(parts)
    return _nonalnum_pattern.sub("", _strip_accents(s))

class APALabelStyle(BaseLabelStyle):
    def format_labels(self, sorted_entries):
        labels = [self.format_label(entry) for entry in sorted_entries]
        count = Counter(labels)
        counted = Counter()
        for label in labels:
            if count[label] == 1:
                yield label
            else:
                yield label + chr(ord('a') + counted[label])
                counted.update([label])

    def format_label(self, entry):
        label = "Anonymous"
        if 'author' in entry.persons:
            label = self.format_author_or_editor_names(entry.persons['author'])
        elif 'editor' in entry.persons:
            label = self.format_author_or_editor_names(entry.persons['editor'])
        elif 'organization' in entry.fields:
            label = entry.fields['organization']
            if label.startswith("The "):
                label = label[4:]

        if 'year' in entry.fields:
            return "{}, {}".format(label, entry.fields['year'])
        else:
            return "{}, n.d.".format(label)

    def format_author_or_editor_names(self, persons):
        if (len(persons) == 1):
            return _strip_nonalnum(persons[0].last_names)
        elif (len(persons) == 2):
            return "{} & {}".format(
                _strip_nonalnum(persons[0].last_names),
                _strip_nonalnum(persons[1].last_names))
        else:
            return "{} et al.".format(
                _strip_nonalnum(persons[0].last_names))

class APAStyle(UnsrtStyle):

    default_label_style = APALabelStyle

register_plugin('pybtex.style.formatting', 'apa', APAStyle)
# # ====================================================================

