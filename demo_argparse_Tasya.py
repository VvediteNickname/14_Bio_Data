import argparse

# ФУНКЦИИ ДЛЯ ВЫПОЛНЕНИЯ КОМАНД

def cmd_fasta(args):
    from sequence_reader import FastaReader
    reader = FastaReader(args.input)
    
    # Если действие не указано (None) — делаем всё
    action = args.action 
    
    if action in (None, "count"):
        print(f"Кол-во последовательностей: {reader.get_sequences_count()}")
    if action in (None, "avg_len"):
        print(f"Средняя длина: {reader.get_average_length()}")
    if action in (None, "print"):
        for name, dna in reader.read():
            print(f">{name}\n{dna}")

def cmd_fastq(args):
    from sequence_reader import FastqReader
    reader = FastqReader(args.input)
    
    action = args.action
    
    if action in (None, "count"):
        print(f"Кол-во последовательностей: {reader.get_sequences_count()}")
    if action in (None, "avg_len"):
        print(f"Средняя длина: {reader.get_average_length()}")
    if action in (None, "plot"):
        reader.plot_quality()
        print("График качества сохранен в 'quality_plot.png'")

def cmd_sam(args):
    from sam_reader import SamReader
    reader = SamReader(args.input)
    
    action = args.action
    
    if action in (None, "count"):
        print(f"Кол-во выравниваний: {reader.count_records()}")
    if action in (None, "stats"):
        print("Статистика по хромосомам:")
        print(reader.chromosome_stats())
    if action == "region":
        print(f"Выравнивания в регионе {args.chrom}:{args.start}-{args.end}:")
        count = 0
        for rec in reader.get_region(args.chrom, args.start, args.end):
            count += 1
        print(f"Найдено: {count}")

def cmd_vcf(args):
    from vcf_reader import VcfReader
    reader = VcfReader(args.input)
    
    action = args.action
    
    if action in (None, "count"):
        print(f"Кол-во вариантов: {reader.count_variants()}")
    if action == "stats":
        print("Статистика по регионам (окно 100000):")
        print(reader.region_count_stats(bin_size=100000))
    if action == "region":
        print(f"Варианты в регионе {args.chrom}:{args.start}-{args.end}:")
        count = 0
        for rec in reader.get_region(args.chrom, args.start, args.end):
            count += 1
        print(f"Найдено: {count}")


# НАСТРОЙКА ПАРСЕРА АРГУМЕНТОВ

def main():
    parser = argparse.ArgumentParser(description="Программа для чтения файлов с биологическими последовательностями.")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Доступные команды")

    # ---------- FASTQ ----------
    p_fastq = subparsers.add_parser("fastq", help="Анализ FASTQ")
    p_fastq.add_argument("-i", "--input", required=True, help="Входной файл")
    p_fastq.add_argument("-a", "--action", choices=["count", "avg_len", "plot"],help="Действие: count (кол-во), avg_len (ср. длина), plot (график)")
    p_fastq.set_defaults(func=cmd_fastq)

    # ---------- FASTA ----------
    p_fasta = subparsers.add_parser("fasta", help="Анализ FASTA")
    p_fasta.add_argument("-i", "--input", required=True, help="Входной файл")
    p_fasta.add_argument("-a", "--action", choices=["count", "avg_len", "print"],help="Действие: count (кол-во), avg_len (ср. длина), print (вывести всё)")
    p_fasta.set_defaults(func=cmd_fasta)

    # ---------- SAM ----------
    p_sam = subparsers.add_parser("sam", help="Анализ SAM")
    p_sam.add_argument("-i", "--input", required=True, help="Входной файл")
    p_sam.add_argument("-a", "--action", choices=["count", "stats", "region"],help="Действие: count (кол-во), stats (статистика), region (поиск в регионе)")
    # Аргументы для поиска в регионе (нужны, только если action == 'region')
    p_sam.add_argument("--chrom", help="Хромосома для поиска (например, chr1)")
    p_sam.add_argument("--start", type=int, help="Начало региона")
    p_sam.add_argument("--end", type=int, help="Конец региона")
    p_sam.set_defaults(func=cmd_sam)

    # ---------- VCF ----------
    p_vcf = subparsers.add_parser("vcf", help="Анализ VCF")
    p_vcf.add_argument("-i", "--input", required=True, help="Входной файл")
    p_vcf.add_argument("-a", "--action", choices=["count", "stats", "region"],help="Действие: count (кол-во), stats (статистика), region (поиск в регионе)")
    p_vcf.add_argument("--chrom", help="Хромосома для поиска")
    p_vcf.add_argument("--start", type=int, help="Начало региона")
    p_vcf.add_argument("--end", type=int, help="Конец региона")
    p_vcf.set_defaults(func=cmd_vcf)

    args = parser.parse_args()
    
    # Проверка: если выбрали 'region', но не указали координаты
    if getattr(args, "action", None) == "region":
        if not (args.chrom and args.start and args.end):
            parser.error("Для действия 'region' необходимо указать --chrom, --start и --end")

    # Вызываем нужную функцию
    args.func(args)

if __name__ == "__main__":
    main()