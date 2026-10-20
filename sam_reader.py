import re
import pandas as pd

from genomic_data_reader import GenomicDataReader


class SamReader(GenomicDataReader):
    """Чтение SAM-файлов."""

    CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")

    def __init__(self, path: str):
        super().__init__(path)
        self._header_loaded = False

    def _load_header(self) -> None:
        """Загружает строки заголовка до первой записи выравнивания."""
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
        """Возвращает идентификаторы записей выбранной группы заголовка."""
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
        """Возвращает первое значение ключа из заголовка, либо None."""
        self._load_header()
        for line in self._meta_lines:
            for field in line.split("\t")[1:]:
                if field.startswith(key + ":"):
                    return field[len(key) + 1:]
        return None

    def _parse_line(self, line: str) -> dict:
        """Разбирает одну строку выравнивания SAM."""
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
        """Последовательно выдаёт записи, не загружая весь файл в память."""
        self._load_header()
        with open(self.path, "r", encoding="utf-8") as file:
            for line in file:
                if not line.strip() or line.startswith("@"):
                    continue
                yield self._parse_line(line)

    def _record_interval(self, record):
        """Возвращает (хромосома, начало, конец) для выравнивания."""
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
        """Возвращает pandas-таблицу с числом выравниваний по хромосомам."""
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
