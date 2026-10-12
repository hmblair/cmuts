# conf.py -- configuration for the documentation site.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

project = "cmuts"
html_title = "cmuts"
author = "Hamish M. Blair"

# The first year is the one in LICENSE. %Y is the year the site is built.
copyright = "2024-%Y, Hamish M. Blair"

extensions = [
    "myst_parser",
    "sphinx_inline_tabs",
]

# ::: fences, {.class} blocks, and $ math.
myst_enable_extensions = ["colon_fence", "attrs_block", "dollarmath"]

# Gives every heading down to this depth an anchor, which the dataset tables link to.
myst_heading_anchors = 3

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# The Makefile renders each version with html_baseurl set to its URL, and passes
# the version list to the sidebar dropdown through the template context.
templates_path = ["_templates"]
html_sidebars = {
    "**": [
        "sidebar/brand.html",
        "sidebar/versions.html",
        "sidebar/search.html",
        "sidebar/scroll-start.html",
        "sidebar/navigation.html",
        "sidebar/scroll-end.html",
    ]
}

# Keeps GitHub Pages from running the site through Jekyll, which drops the
# _static and _sources directories.
html_extra_path = [".nojekyll"]

# The Makefile sets source_branch to the ref each version is rendered from.
html_theme_options = {
    "source_repository": "https://github.com/DasLab/cmuts",
    "source_branch": "main",
    "source_directory": "docs/",
}
