'''
from fasta_reader import FastaReader
from fastq_reader import FastqReader
from sam_reader import SamReader
'''
from vcf_reader import VcfReader


'''
    Общая красивенькая часть
'''
print("=" * 33)
print("Выберите формат тестируемого файла (fasta, fastq, sam, vcf): ")
file_format = input()




'''
    Часть для vcf-файлов.
'''
if file_format in "vcf":
    print("Выберите файл для тестирования: sample.vcf, 88_1_9001000-10001000_fullgenome.vcf")
    filename = input()
    print("=" * 33)
    
    r = VcfReader(filename)
    r._load_header()

    print("Мета-строки:")
    for m in r.get_meta_lines():
        print("  ", m)

    print("\nINFO ID:", " ".join(r.get_group_ids("INFO")))
    print("FORMAT ID:", " ".join(r.get_group_ids("FORMAT")))
    print("FILTER ID:", " ".join(r.get_group_ids("FILTER")))
    print("contig ID:", " ".join(r.get_group_ids("contig")))

    print("\nfileformat:", r.get_meta_value("fileformat"))
    print("source:", r.get_meta_value("source"))
    print("contig (простой ключ):", r.get_meta_value("contig"))

    print("\nКолонки:", " ".join(r.get_columns()))

    print("Кол-во вариантов:", r.count_variants())

    print("\nСтатистика по регионам (окно 100000):")
    stats = r.region_count_stats(bin_size=100000)
    print(stats)

    n = sum(1 for _ in r.get_region("1", 9000001, 9100000))
    print(f"\nВариантов в 1:9000001-9100000 = {n}")