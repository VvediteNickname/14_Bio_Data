from abc import ABC, abstractmethod
from reader import Reader   


class GenomicDataReader(Reader):  
    """Базовый класс: общие поля, общий заголовок, общее количество записей,
    общий поиск по региону. Формат-специфика — в наследниках."""

    def __init__(self, path: str):
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
        """
        Реализация абстрактного метода из Reader.
        Перенаправляет вызов на __iter__, который реализован в наследниках.
        """
        return self.__iter__()

    # ---------- ОБЩЕЕ: одинаково у VCF и SAM ----------

    def get_meta_lines(self) -> list[str]:
        """Все строки заголовка как есть (##... / @...)."""
        self._load_header()
        return list(self._meta_lines)

    def get_header_line(self) -> str:
        """Строка с именами колонок (#CHROM ... у VCF; у SAM может быть пустой)."""
        self._load_header()
        return self._header_line

    def get_columns(self) -> list[str]:
        """Колонки из строки #CHROM ... (у SAM может быть пустой список)."""
        self._load_header()
        return list(self._columns)

    def count_records(self) -> int:
        """Общее число записей (вариантов / выравниваний). Память O(1)."""
        return sum(1 for _ in self)

    def get_region(self, chrom, start, end, overlap=True):
        """Аналог bedtools intersect: записи, пересекающиеся с [start, end]
        (1-based, включительно). Как именно считать интервал записи —
        решает наследник в _record_interval()."""
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
        """Прочитать заголовок файла и заполнить _meta_lines / _header_line / _columns."""
        ...

    @abstractmethod
    def __iter__(self):
        """Ленивая итерация по записям формата (dict)."""
        ...

    @abstractmethod
    def _record_interval(self, record):
        """Вернуть (chrom, start, end) для одной записи.
        Координаты 1-based, включительно."""
        ...

    @abstractmethod
    def get_group_ids(self, group: str):
        """ID внутри группы заголовка.
        VCF: '##INFO=<ID=DP,...>'  -> 'DP'
        SAM: '@SQ\\tSN:chr1\\t...' -> 'chr1'"""
        ...

    @abstractmethod
    def get_meta_value(self, key: str):
        """Значение простого ключа заголовка.
        VCF: '##source=myCaller' -> 'myCaller'
        SAM: '@HD\\tVN:1.6'      -> 'VN:1.6'"""
        ...