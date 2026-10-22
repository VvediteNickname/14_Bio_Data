from abc import ABC, abstractmethod
from reader import Reader   


class GenomicDataReader(Reader):  
    """Базовый класс для чтения геномных форматов (SAM, VCF).

    Содержит общие поля, общую работу с заголовком, подсчёт числа записей
    и поиск по региону. Формат-специфичная логика реализуется в наследниках.

    Attributes:
        path (str): Путь к файлу.
        _meta_lines (list[str]): Строки заголовка как есть.
        _header_line (str): Строка с именами колонок.
        _columns (list[str]): Имена колонок из строки заголовка.
        _header_loaded (bool): Флаг, указывающий, загружен ли заголовок.
    """

    def __init__(self, path: str):
        """Инициализирует читатель геномных данных.

        Args:
            path (str): Путь к файлу.
        """
        # Передаем путь в родительский класс Reader
        super().__init__(path)  
        
        # Оставляем self.path, чтобы не переписывать весь код в SamReader и VcfReader
        self.path = path  
        
        self._meta_lines: list[str] = []
        self._header_line: str = ""
        self._columns: list[str] = []
        self._header_loaded = False  # ← Полезно вынести сюда, чтобы не дублировать в SAM и VCF

    # ← ДОБАВЛЕНО: мостик между контрактом Reader (read) и итерацией (__iter__)
    def read(self):
        """Реализация абстрактного метода из Reader.

        Перенаправляет вызов на :meth:`__iter__`, который реализован
        в наследниках.

        Returns:
            Iterator: Итератор по записям формата.
        """
        return self.__iter__()

    # ---------- ОБЩЕЕ: одинаково у VCF и SAM ----------

    def get_meta_lines(self) -> list[str]:
        """Возвращает все строки заголовка как есть.

        Returns:
            list[str]: Копия списка строк заголовка (``##...`` / ``@...``).
        """
        self._load_header()
        return list(self._meta_lines)

    def get_header_line(self) -> str:
        """Возвращает строку с именами колонок.

        Returns:
            str: Строка заголовка (``#CHROM ...`` у VCF; у SAM может быть
            пустой).
        """
        self._load_header()
        return self._header_line

    def get_columns(self) -> list[str]:
        """Возвращает имена колонок из строки заголовка.

        Returns:
            list[str]: Копия списка колонок (у SAM может быть пустой список).
        """
        self._load_header()
        return list(self._columns)

    def count_records(self) -> int:
        """Возвращает общее число записей.

        Подсчитывает варианты (VCF) или выравнивания (SAM), не загружая
        весь файл в память (память O(1)).

        Returns:
            int: Число записей.
        """
        return sum(1 for _ in self)

    def get_region(self, chrom, start, end, overlap=True):
        """Итерирует по записям, пересекающимся с заданным регионом.

        Аналог ``bedtools intersect``. Координаты 1-based, границы включены.
        Способ вычисления интервала записи определяется наследником
        в :meth:`_record_interval`.

        Args:
            chrom (str): Имя хромосомы.
            start (int): Начало региона (1-based, включительно).
            end (int): Конец региона (1-based, включительно).
            overlap (bool, optional): Если ``True``, возвращаются записи,
                пересекающиеся с регионом; если ``False`` — полностью
                попадающие внутрь региона. По умолчанию ``True``.

        Yields:
            dict: Очередная запись, удовлетворяющая условию.

        Raises:
            ValueError: Если ``start > end``.
        """
        if start > end:
            raise ValueError("start > end")
        for rec in self:
            r_chrom, r_start, r_end = self._record_interval(rec)
            if r_chrom != chrom:
                continue
            if overlap:
                if r_start <= end and r_end >= start:
                    yield rec
            else:
                if r_start >= start and r_end <= end:
                    yield rec

    # ---------- АБСТРАКТНОЕ: у каждого формата своё ----------

    @abstractmethod
    def _load_header(self) -> None:
        """Читает заголовок файла.

        Заполняет ``_meta_lines``, ``_header_line`` и ``_columns``.

        Returns:
            None
        """
        ...

    @abstractmethod
    def __iter__(self):
        """Ленивая итерация по записям формата.

        Yields:
            dict: Очередная разобранная запись.
        """
        ...

    @abstractmethod
    def _record_interval(self, record):
        """Возвращает интервал одной записи на референсе.

        Args:
            record (dict): Разобранная запись формата.

        Returns:
            tuple[str, int, int]: Кортеж ``(chrom, start, end)``.
            Координаты 1-based, границы включены.
        """
        ...

    @abstractmethod
    def get_group_ids(self, group: str):
        """Возвращает ID внутри указанной группы заголовка.

        Примеры:
            VCF: ``##INFO=<ID=DP,...>`` -> ``'DP'``
            SAM: ``@SQ\\tSN:chr1\\t...`` -> ``'chr1'``

        Args:
            group (str): Имя группы заголовка.

        Returns:
            list[str]: Список идентификаторов группы.
        """
        ...

    @abstractmethod
    def get_meta_value(self, key: str):
        """Возвращает значение простого ключа заголовка.

        Примеры:
            VCF: ``##source=myCaller`` -> ``'myCaller'``
            SAM: ``@HD\\tVN:1.6`` -> ``'VN:1.6'``

        Args:
            key (str): Имя ключа метаданных.

        Returns:
            str | None: Значение ключа или ``None``, если ключ не найден.
        """
        ...
