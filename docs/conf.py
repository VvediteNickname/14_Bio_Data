project = "Программа для работы с основными биологическими форматами"
author = "Команда проекта"
copyright = "2026, Команда проекта"
release = "1.0"

extensions = []

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
html_title = "Программа для работы с основными биологическими форматами"
html_short_title = "Биологические форматы"
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
