Документация Sphinx для проекта «Ква-ква-ква»

В архиве находится готовая папка docs. Скопируйте её в корневой каталог локального проекта, заменив предыдущую папку docs.

Linux / Ubuntu:

1. Откройте терминал в корне проекта.
2. Активируйте виртуальное окружение проекта. Если его нет:
   python3 -m venv .venv
   source .venv/bin/activate
3. Установите библиотеки проекта:
   python -m pip install pandas matplotlib
4. Установите Sphinx и тему:
   python -m pip install -r docs/requirements.txt
5. Соберите HTML-документацию:
   python -m sphinx -b html docs docs/_build/html
6. Откройте сайт:
   xdg-open docs/_build/html/index.html

Внешний вид оформлен в пастельной мятной гамме, с округлёнными блоками и боковыми пиксельными лианами. Навигация формируется Sphinx из разделов index.rst.
