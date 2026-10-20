from sequence_reader import FastaReader, FastqReader
from vcf_reader import VcfReader
from sam_reader import SamReader
from sequence_reader import FastaReader, FastqReader

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

'''
    Часть для sam-файлов.
'''

if file_format == "sam":

    print("Выберите файл для тестирования: sample.sam, example.sam")
    filename = input().strip()

    print("=" * 33)

    r = SamReader(filename)

    r._load_header()

    print("Мета-строки:")

    for m in r.get_meta_lines():
        print(" ", m)

    print("\nPG ID:", " ".join(r.get_group_ids("PG")))
    print("RG ID:", " ".join(r.get_group_ids("RG")))
    print("SQ ID:", " ".join(r.get_group_ids("SQ")))

    print("\nКол-во выравниваний:", r.count_records())

    print("\nСтатистика по хромосомам:")

    stats = r.chromosome_stats()
    print(stats)

    print("\nВыравнивания в регионе chr1:100-500:")

    n = sum(1 for _ in r.get_region("chr1", 100, 500))

    print(f"Кол-во выравниваний в chr1:100-500 = {n}")

    # ГАЙД: как запустить SAM - файл:
    #python3 -m venv .venv
    #source .venv/bin/activate
    #python -m pip install pandas

<<<<<<< HEAD


'''
    ЧАСТЬ ДЛЯ FASTA
'''
if file_format == "fasta":
    print("Выберите файл для тестирования:")
    filename = input()
    
    # создаем робота для FASTA
    r = FastaReader(filename)
    
    # считаем количество и выводим
    count = r.get_sequences_count()
    print("Количество последовательностей в файле:")
    print(count)
    
    # считаем среднюю длину и выводим
    avg_len = r.get_average_length()
    print("Средняя длина последовательностей:")
    print(avg_len)


'''
    ЧАСТЬ ДЛЯ FASTQ
'''
if file_format == "fastq":
    print("Выберите файл для тестирования:")
    filename = input()
    
    # создаем робота для FASTQ
    r = FastqReader(filename)
    
    # считаем количество и выводим
    count = r.get_sequences_count()
    print("Количество последовательностей в файле:")
    print(count)
    
    # считаем среднюю длину и выводим
    avg_len = r.get_average_length()
    print("Средняя длина последовательностей:")
    print(avg_len)
    
    # запускаем создание картинки с графиком
    print("Строим и сохраняем график качества...")
    r.plot_quality()

=======
'''
    Часть для fasta-файлов.
'''

if file_format == "fasta":

    print("Выберите файл для тестирования: sample.fasta")
    filename = input().strip()
    print("=" * 33)

    r = FastaReader(filename)

    print("Общая статистика по файлу:")
    print("Кол-во последовательностей:", r.get_sequences_count())
    print("Средняя длина последовательностей:", r.get_average_length())

    print("\nВывод первых трех последовательностей:")
    count = 0
    for name, dna in r.read():
        # Выводим только начало последовательности, чтобы не засорять консоль
        short_dna = dna[:50] + "..." if len(dna) > 50 else dna
        print(f"  >{name}")
        print(f"  {short_dna}")
        count += 1
        if count >= 3:
            print("  ... (остальные последовательности скрыты для экономии места)")
            break


'''
    Часть для fastq-файлов.
'''

if file_format == "fastq":

    print("Выберите файл для тестирования: sample.fastq")
    filename = input().strip()
    print("=" * 33)

    r = FastqReader(filename)

    print("Общая статистика по файлу:")
    print("Кол-во последовательностей:", r.get_sequences_count())
    print("Средняя длина последовательностей:", r.get_average_length())

    print("\nВывод первых трех последовательностей (с качеством):")
    count = 0
    for name, dna, quality in r.read():
        short_dna = dna[:50] + "..." if len(dna) > 50 else dna
        short_qual = quality[:50] + "..." if len(quality) > 50 else quality
        print(f"  @{name}")
        print(f"  {short_dna}")
        print(f"  +")
        print(f"  {short_qual}")
        count += 1
        if count >= 3:
            print("  ... (остальные последовательности скрыты для экономии места)")
            break

    print("\nПостроение графика качества (Per base sequence quality)...")
    r.plot_quality()
>>>>>>> 11f9ca0cc1283072b4d8352dfc90ff163215f229
