# Configuration file for the Sphinx documentation builder.

import os
import sys
import subprocess

from m2r import MdInclude
from recommonmark.transform import AutoStructify

# check if we are running on readthedocs.org
on_readthedocs = os.environ.get("READTHEDOCS", None) == "True"

# -- Project information ------------------------------------------------------

project = "Acts"
author = "The Acts authors"
copyright = "2014–2020 CERN for the benefit of the Acts project"
# version = '@PROJECT_VERSION@'
# release = '@PROJECT_VERSION@'

# -- General ------------------------------------------------------------------

extensions = [
    "breathe",
    "recommonmark",
    "sphinx.ext.mathjax",
    "sphinx_markdown_tables",
]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
# ensure the in-source build directory is ignored
exclude_patterns = [
    "_build",
]
# cpp as default language
primary_domain = "cpp"
highlight_language = "cpp"
smartquotes = True
numfig = True

# -- Options for HTML output --------------------------------------------------

# ensure we use the RTD them when building locally
if not on_readthedocs:
    import sphinx_rtd_theme

    html_theme = "sphinx_rtd_theme"
    html_theme_path = [
        sphinx_rtd_theme.get_html_theme_path(),
    ]

html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 2,
    "prev_next_buttons_location": None,  # no prev/next links
    "style_external_links": True,
}
html_logo = "figures/acts_logo_white.svg"
html_static_path = [
    "_static",
]
html_css_files = [
    "custom.css",
]
html_copy_source = False
html_show_sourcelink = False
html_show_sphinx = False

# -- Doxygen integration with Breathe -----------------------------------------

breathe_projects = {
    "Acts": "_build/doxygen-xml",
}
breathe_default_project = "Acts"
breathe_domain_by_extension = {
    "cpp": "cpp",
    "hpp": "cpp",
    "ipp": "cpp",
}
breathe_default_members = (
    "members",
    "undoc-members",
)

# -- Automatic API documentation ---------------------------------------------

env = os.environ.copy()
env["DOXYGEN_WARN_AS_ERROR"] = "NO"
cwd = os.path.dirname(__file__)

if on_readthedocs or tags.has("run_doxygen"):
    # if we are running on RTD Doxygen must be run as part of the build
    print("Executing doxygen in", cwd)
    print(
        "Doxygen version:",
        subprocess.check_output(["doxygen", "--version"], encoding="utf-8"),
    )
    sys.stdout.flush()
    subprocess.check_call(
        ["doxygen", "Doxyfile"], stdout=subprocess.PIPE, cwd=cwd, env=env
    )

if on_readthedocs or tags.has("run_apidoc"):
    print("Executing breathe apidoc in", cwd)
    subprocess.check_call(
        [sys.executable, "-m", "breathe.apidoc", "_build/doxygen-xml", "-o", "api"],
        stdout=subprocess.PIPE,
        cwd=cwd,
        env=env,
    )

# -- Markdown bridge setup hook (must come last, not sure why) ----------------


def setup(app):
    app.add_config_value(
        "recommonmark_config",
        {
            "enable_math": True,
            "enable_inline_math": True,
        },
        True,
    )
    app.add_transform(AutoStructify)

    app.add_config_value("no_underscore_emphasis", False, "env")
    app.add_config_value("m2r_parse_relative_links", False, "env")
    app.add_config_value("m2r_anonymous_references", False, "env")
    app.add_config_value("m2r_disable_inline_math", False, "env")
    app.add_directive("mdinclude", MdInclude)
