import re
import pandas as pd
from genomic_data_reader import GenomicDataReader


class VcfReader(GenomicDataReader):
    """Чтение VCF-файлов.

    Класс обеспечивает разбор заголовка VCF, извлечение метаданных,
    последовательную итерацию по вариантам, подсчёт статистики по регионам
    и совместимость с интерфейсом :class:`GenomicDataReader`.

    Attributes:
        path (str): Путь к VCF-файлу (наследуется от GenomicDataReader).
        _header_loaded (bool): Флаг, указывающий, загружен ли заголовок.
        _n_header_lines (int): Число строк заголовка (``##`` и ``#CHROM``).
        _meta_lines (list[str]): Строки метаданных ``##``.
        _header_line (str): Строка ``#CHROM ...``.
        _columns (list[str]): Имена колонок из строки ``#CHROM``.
    """

    def __init__(self, path: str):
        """Инициализирует читатель VCF-файла.

        Args:
            path (str): Путь к VCF-файлу.
        """
        super().__init__(path)
        self._header_loaded = False
        self._n_header_lines = 0

    # ---------- заголовок ----------
    def _load_header(self) -> None:
        """Загружает заголовок VCF-файла.

        Читает строки метаданных ``##`` в ``self._meta_lines`` и строку
        ``#CHROM`` в ``self._header_line`` / ``self._columns``, подсчитывая
        общее число строк заголовка в ``self._n_header_lines``.
        Повторный вызов не перечитывает файл.

        Returns:
            None

        Raises:
            ValueError: Если данные встречаются до строки ``#CHROM`` или
                строка ``#CHROM`` отсутствует (файл пуст или не является VCF).
        """
        if self._header_loaded:
            return
        self._n_header_lines = 0
        with open(self.path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("##"):
                    self._meta_lines.append(line)
                    self._n_header_lines += 1
                elif line.startswith("#"):
                    self._header_line = line
                    self._columns = line.lstrip("#").split("\t")
                    self._n_header_lines += 1
                    break
                else:
                    raise ValueError("Данные до заголовка #CHROM — это не VCF")
            else:
                # проверка [в файле вообще нет строки #CHROM — это не VCF (или файл пуст)]
                raise ValueError(
                    f"VCF: в файле {self.path} нет строки #CHROM — это не VCF или файл пуст"
                )
        self._header_loaded = True

    # ---------- VCF-специфичные парсеры заголовка ----------
    def get_group_ids(self, group):
        """Возвращает список ID для указанной группы метаданных VCF.

        Ищет строки вида ``##<group>=<ID=...>`` и извлекает значения ``ID``.

        Args:
            group (str): Имя группы метаданных (например, ``"INFO"``,
                ``"FILTER"``, ``"FORMAT"``, ``"contig"``).

        Returns:
            list[str]: Список значений ``ID`` в порядке следования строк.

        Raises:
            ValueError: Если ``group`` не является непустой строкой.
        """

        # проверка [group — непустая строка]
        if not isinstance(group, str) or not group:
            raise ValueError("group должен быть непустой строкой")
        
        self._load_header()
        pattern = re.compile(rf"^##{group}=<ID=([^,>]+)")
        result = []
        for line in self._meta_lines:
            m = pattern.match(line)
            if m:
                result.append(m.group(1))
        return result

    def get_meta_value(self, key):
        """Возвращает первое значение метаданных по ключу.

        Ищет строку вида ``##<key>=<value>`` и возвращает часть после ``=``.

        Args:
            key (str): Имя ключа метаданных (например, ``"fileformat"``,
                ``"reference"``).

        Returns:
            str | None: Значение ключа или ``None``, если ключ не найден.

        Raises:
            ValueError: Если ``key`` не является непустой строкой.
        """

        # проверка [key — непустая строка]
        if not isinstance(key, str) or not key:
            raise ValueError("key должен быть непустой строкой")
        
        self._load_header()
        prefix = f"##{key}="
        for line in self._meta_lines:
            if line.startswith(prefix):
                return line[len(prefix):]
        return None

    # ---------- парсинг одной строки данных ----------
    def _parse_line(self, line):
        """Разбирает одну строку данных VCF.

        Извлекает обязательные колонки, разбирает поле ``INFO`` в словарь
        и приводит ``QUAL`` к ``float`` (или ``None`` для ``"."``). Поля
        ``FORMAT`` и образцы сохраняются как есть.

        Args:
            line (str): Строка данных VCF.

        Returns:
            dict: Словарь с ключами ``chrom``, ``pos``, ``id``, ``ref``,
            ``alt``, ``qual``, ``filt``, ``info``, ``format``, ``samples``.

        Raises:
            ValueError: Если в строке меньше 8 колонок, ``REF`` пустой или
                равен ``"."``, либо ``POS`` не является положительным целым.
        """
        fields = line.rstrip("\n").split("\t")

        # проверка [в строке минимум 8 колонок — иначе нечего распаковывать]
        if len(fields) < 8:
            raise ValueError(
                f"VCF: ожидалось минимум 8 колонок, получено {len(fields)}: "
                f"{line[:80]!r}"
            )

        chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]

        # проверка [REF не пустой и не '.', иначе _record_interval даст end < pos]
        if not ref or ref == ".":
            raise ValueError(f"VCF: REF пустой в позиции {chrom}:{pos}")

        # проверка [POS — целое положительное число]
        pos_int = int(pos)  # ValueError, если не число
        if pos_int <= 0:
            raise ValueError(f"VCF: POS должен быть > 0, получено {pos_int} в {chrom}")

        info_dict = {}
        if info and info != ".":
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_dict[k] = v
                else:
                    info_dict[kv] = True

        try:
            qual_f = None if qual == "." else float(qual)
        except ValueError:
            qual_f = None

        fmt = fields[8] if len(fields) > 8 else None
        samples = fields[9:] if len(fields) > 9 else []

        return {
            "chrom": chrom, "pos": int(pos), "id": vid,
            "ref": ref, "alt": alt, "qual": qual_f,
            "filt": filt, "info": info_dict,
            "format": fmt, "samples": samples,
        }

    # ---------- итерация ----------
    def __iter__(self):
        """Последовательно выдаёт варианты VCF.

        Пропускает строки заголовка и пустые строки. Каждая строка данных
        разбирается методом :meth:`_parse_line`.

        Yields:
            dict: Очередной разобранный вариант VCF.
        """
        self._load_header()
        with open(self.path, "r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                if i < self._n_header_lines:
                    continue
                if not line.strip() or line.startswith("#"):
                    continue
                yield self._parse_line(line)

    # ---------- интервал варианта (для get_region из родителя) ----------
    def _record_interval(self, v):
        """Возвращает интервал варианта на референсе.

        Args:
            v (dict): Разобранный вариант VCF.

        Returns:
            tuple[str, int, int]: Кортеж ``(хромосома, начало, конец)``,
            где конец вычисляется как ``pos + len(ref) - 1``.
        """
        return (v["chrom"], v["pos"], v["pos"] + len(v["ref"]) - 1)

    # ---------- «человеческие» обёртки ----------
    def count_variants(self):
        """Возвращает общее число вариантов в файле.

        Является удобной обёрткой над :meth:`count_records`.

        Returns:
            int: Число вариантов.
        """
        return self.count_records()

    # ---------- статистика (VCF-специфична) ----------
    def region_count_stats(self, bin_size=100000, chrom=None):
        """Считает число вариантов по бинам заданного размера.

        Разбивает координаты на бины шириной ``bin_size`` (1-based) и
        подсчитывает число вариантов в каждом бине. При необходимости
        ограничивает подсчёт одной хромосомой.

        Args:
            bin_size (int, optional): Размер бина в парах оснований.
                Должен быть положительным целым. По умолчанию ``100000``.
            chrom (str | None, optional): Если задано, учитываются только
                варианты указанной хромосомы. По умолчанию ``None``.

        Returns:
            pd.DataFrame: Таблица с колонками ``chrom``, ``start``, ``end``,
            ``count``. Пустая таблица с теми же колонками, если вариантов нет.

        Raises:
            ValueError: Если ``bin_size`` не является положительным целым.
        """

        # проверка [bin_size — положительное целое, иначе ZeroDivisionError]
        if not isinstance(bin_size, int) or bin_size <= 0:
            raise ValueError(
                f"bin_size должен быть положительным целым, получено {bin_size!r}"
            )

        rows = []
        for v in self:
            if chrom is not None and v["chrom"] != chrom:
                continue
            bin_start = (v["pos"] - 1) // bin_size * bin_size + 1
            rows.append((v["chrom"], bin_start))

        if not rows:
            return pd.DataFrame(columns=["chrom", "start", "end", "count"])

        df = pd.DataFrame(rows, columns=["chrom", "start"])
        stats = (df.groupby(["chrom", "start"])
                   .size()
                   .reset_index(name="count"))
        stats["end"] = stats["start"] + bin_size - 1
        return stats[["chrom", "start", "end", "count"]]
        
