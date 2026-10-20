project = "Программа для работы с основными биологическими форматами"
author = "Ква-ква-ква"
copyright = "2026, Ква-ква-ква"
release = "1.0"

extensions = ["sphinx.ext.graphviz"]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
html_title = "Ква-ква-ква — документация проекта"
html_short_title = "Ква-ква-ква"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 3,
    "includehidden": True,
    "titles_only": False,
    "style_nav_header_background": "#A9DFC0",
}

pygments_style = "friendly"
