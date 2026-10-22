import re
import pandas as pd

from genomic_data_reader import GenomicDataReader


class SamReader(GenomicDataReader):
    """Чтение SAM-файлов.

    Класс предоставляет последовательный доступ к записям выравнивания SAM,
    разбор заголовка, извлечение метаданных и подсчёт статистики по хромосомам.

    Attributes:
        CIGAR_PATTERN (re.Pattern): Регулярное выражение для разбора CIGAR-строк.
        path (str): Путь к SAM-файлу (наследуется от GenomicDataReader).
        _header_loaded (bool): Флаг, указывающий, загружен ли заголовок.
        _meta_lines (list[str]): Строки заголовка SAM-файла.
    """

    CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")

    def __init__(self, path: str):
        """Инициализирует читатель SAM-файла.

        Args:
            path (str): Путь к SAM-файлу.
        """
        super().__init__(path)
        self._header_loaded = False

    def _load_header(self) -> None:
        """Загружает строки заголовка до первой записи выравнивания.

        Читает файл построчно и сохраняет все строки, начинающиеся с ``@``,
        в ``self._meta_lines``. Останавливается на первой строке, не являющейся
        частью заголовка. Повторный вызов не перечитывает файл.

        Returns:
            None
        """
        if self._header_loaded:
            return

        self._meta_lines = []

        with open(self.path, "r", encoding="utf-8") as file:
            for line in file:
                line = line.rstrip("\r\n")
                if line.startswith("@"):
                    self._meta_lines.append(line)
                else:
                    break

        self._header_loaded = True

    def get_group_ids(self, group: str) -> list[str]:
        """Возвращает идентификаторы записей выбранной группы заголовка.

        Для групп ``SQ`` и ``HD`` ищет поле ``SN`` или ``ID`` соответственно.
        Для группы ``CO`` возвращает текст комментария. Для прочих групп
        используется ключ ``ID``.

        Args:
            group (str): Имя группы заголовка (например, ``"SQ"``, ``"HD"``,
                ``"CO"``). Ведущий символ ``@`` игнорируется.

        Returns:
            list[str]: Список идентификаторов (или текстов комментариев для
            группы ``CO``) в порядке следования строк заголовка.
        """
        self._load_header()
        group = group.lstrip("@")
        result = []
        id_key = "SN" if group == "SQ" else "ID"

        for line in self._meta_lines:
            fields = line.split("\t")
            if fields[0] != f"@{group}":
                continue

            if group == "CO":
                result.append("\t".join(fields[1:]))
                continue

            for field in fields[1:]:
                if field.startswith(id_key + ":"):
                    result.append(field[len(id_key) + 1:])
                    break

        return result

    def get_meta_value(self, key: str):
        """Возвращает первое значение ключа из заголовка, либо None.

        Просматривает все поля всех строк заголовка и возвращает значение
        первого поля, начинающегося с ``key + ":"``.

        Args:
            key (str): Имя ключа метаданных (например, ``"VN"``, ``"SO"``).

        Returns:
            str | None: Значение ключа без префикса или ``None``, если ключ
            не найден.
        """
        self._load_header()
        for line in self._meta_lines:
            for field in line.split("\t")[1:]:
                if field.startswith(key + ":"):
                    return field[len(key) + 1:]
        return None

    def _parse_line(self, line: str) -> dict:
        """Разбирает одну строку выравнивания SAM.

        Преобразует обязательные поля записи SAM в словарь, приводя числовые
        поля к ``int``. Дополнительные поля сохраняются в ключе ``tags``
        без преобразования.

        Args:
            line (str): Строка SAM-файла с записью выравнивания.

        Returns:
            dict: Словарь с ключами ``name``, ``flag``, ``chrom``, ``pos``,
            ``mapq``, ``cigar``, ``rnext``, ``pnext``, ``tlen``, ``seq``,
            ``qual`` и ``tags``.

        Raises:
            ValueError: Если в строке меньше 11 полей или числовое поле
                некорректно.
        """
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise ValueError("В записи SAM должно быть не менее 11 полей")

        names = (
            "name", "flag", "chrom", "pos", "mapq",
            "cigar", "rnext", "pnext", "tlen", "seq", "qual"
        )
        record = dict(zip(names, fields[:11]))

        try:
            record["flag"] = int(record["flag"])
            record["pos"] = int(record["pos"])
            record["mapq"] = int(record["mapq"])
            record["pnext"] = (
                None if record["pnext"] == "*"
                else int(record["pnext"])
            )
            record["tlen"] = int(record["tlen"])
        except ValueError as exc:
            raise ValueError("Некорректное числовое поле в записи SAM") from exc

        # Дополнительные поля сохраняются без преобразования.
        record["tags"] = fields[11:]
        return record

    def __iter__(self):
        """Последовательно выдаёт записи, не загружая весь файл в память.

        Пропускает пустые строки и строки заголовка. Каждая запись
        разбирается методом :meth:`_parse_line`.

        Yields:
            dict: Очередная разобранная запись выравнивания SAM.
        """
        self._load_header()
        with open(self.path, "r", encoding="utf-8") as file:
            for line in file:
                if not line.strip() or line.startswith("@"):
                    continue
                yield self._parse_line(line)

    def _record_interval(self, record):
        """Возвращает интервал выравнивания на референсе.

        Вычисляет длину участка на референсе по CIGAR-строке, учитывая
        операции ``M``, ``D``, ``N``, ``=`` и ``X``.

        Args:
            record (dict): Разобранная запись выравнивания SAM.

        Returns:
            tuple[str, int, int]: Кортеж ``(хромосома, начало, конец)``.
            Для несопоставленных записей (``chrom == "*"``, ``pos == 0``
            или ``cigar == "*"``) возвращается ``(chrom, start, start)``.

        Raises:
            ValueError: Если CIGAR-строка некорректна или не задаёт длину
                на референсе.
        """
        chrom = record["chrom"]
        start = record["pos"]
        cigar = record["cigar"]

        if chrom == "*" or start == 0 or cigar == "*":
            return chrom, start, start

        operations = self.CIGAR_PATTERN.findall(cigar)
        reconstructed = "".join(length + operation for length, operation in operations)
        if not operations or reconstructed != cigar:
            raise ValueError(f"Некорректный CIGAR: {cigar}")

        # Эти операции занимают позиции на референсе.
        reference_length = sum(
            int(length)
            for length, operation in operations
            if operation in "MDN=X"
        )
        if reference_length == 0:
            raise ValueError(f"CIGAR не задаёт длину на референсе: {cigar}")

        # В SAM координаты 1-based, обе границы включены.
        return chrom, start, start + reference_length - 1

    def chromosome_stats(self) -> pd.DataFrame:
        """Возвращает pandas-таблицу с числом выравниваний по хромосомам.

        Пропускает несопоставленные записи (``chrom == "*"`` или ``pos == 0``).
        Результат отсортирован по имени хромосомы.

        Returns:
            pd.DataFrame: Таблица с колонками ``chrom`` и ``alignments``,
            содержащая число выравниваний для каждой хромосомы.
        """
        counts = {}
        for record in self:
            chrom = record["chrom"]
            if chrom == "*" or record["pos"] == 0:
                continue
            counts[chrom] = counts.get(chrom, 0) + 1

        return pd.DataFrame(
            [{"chrom": chrom, "alignments": count}
             for chrom, count in sorted(counts.items())],
            columns=["chrom", "alignments"],
        )
