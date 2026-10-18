import re
import pandas as pd
from genomic_data_reader import GenomicDataReader


class VcfReader(GenomicDataReader):

    def __init__(self, path: str):
        super().__init__(path)
        self._header_loaded = False
        self._n_header_lines = 0

    # ---------- заголовок ----------
    def _load_header(self) -> None:
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
        self._header_loaded = True

    # ---------- VCF-специфичные парсеры заголовка ----------
    def get_group_ids(self, group):
        self._load_header()
        pattern = re.compile(rf"^##{group}=<ID=([^,>]+)")
        result = []
        for line in self._meta_lines:
            m = pattern.match(line)
            if m:
                result.append(m.group(1))
        return result

    def get_meta_value(self, key):
        self._load_header()
        prefix = f"##{key}="
        for line in self._meta_lines:
            if line.startswith(prefix):
                return line[len(prefix):]
        return None

    # ---------- парсинг одной строки данных ----------
    def _parse_line(self, line):
        fields = line.rstrip("\n").split("\t")
        chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]

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
        return (v["chrom"], v["pos"], v["pos"] + len(v["ref"]) - 1)

    # ---------- «человеческие» обёртки ----------
    def count_variants(self):
        return self.count_records()

    # ---------- статистика (VCF-специфична) ----------
    def region_count_stats(self, bin_size=100000, chrom=None):
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