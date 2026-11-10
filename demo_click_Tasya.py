import click

# ГРУППА КОМАНД (Главная точка входа)
@click.group()
def main():
    """Программа для чтения файлов с биологическими последовательностями."""
    pass


# КОМАНДА: FASTA

@main.command()
@click.option("-i", "--input", "input_file", required=True, help="Входной файл")
@click.option("-a", "--action", type=click.Choice(["count", "avg_len", "print"]), default=None,
              help="Действие: count (кол-во), avg_len (ср. длина), print (вывести всё)")
def fasta(input_file, action):
    """Анализ FASTA"""
    from sequence_reader import FastaReader
    reader = FastaReader(input_file)
    
    if action in (None, "count"):
        click.echo(f"Кол-во последовательностей: {reader.get_sequences_count()}")
    if action in (None, "avg_len"):
        click.echo(f"Средняя длина: {reader.get_average_length()}")
    if action in (None, "print"):
        for name, dna in reader.read():
            click.echo(f">{name}\n{dna}")


# КОМАНДА: FASTQ

@main.command()
@click.option("-i", "--input", "input_file", required=True, help="Входной файл")
@click.option("-a", "--action", type=click.Choice(["count", "avg_len", "plot"]), default=None,
              help="Действие: count (кол-во), avg_len (ср. длина), plot (график)")
def fastq(input_file, action):
    """Анализ FASTQ"""
    from sequence_reader import FastqReader
    reader = FastqReader(input_file)
    
    if action in (None, "count"):
        click.echo(f"Кол-во последовательностей: {reader.get_sequences_count()}")
    if action in (None, "avg_len"):
        click.echo(f"Средняя длина: {reader.get_average_length()}")
    if action in (None, "plot"):
        reader.plot_quality()
        click.echo("График качества сохранен в 'quality_plot.png'")


# КОМАНДА: SAM

@main.command()
@click.option("-i", "--input", "input_file", required=True, help="Входной файл")
@click.option("-a", "--action", type=click.Choice(["count", "stats", "region"]), default=None,
              help="Действие: count (кол-во), stats (статистика), region (поиск в регионе)")
@click.option("--chrom", help="Хромосома для поиска (например, chr1)")
@click.option("--start", type=int, help="Начало региона")
@click.option("--end", type=int, help="Конец региона")
def sam(input_file, action, chrom, start, end):
    """Анализ SAM"""
    # Проверка: если выбрали 'region', но не указали координаты
    if action == "region" and not (chrom and start and end):
        raise click.UsageError("Для действия 'region' необходимо указать --chrom, --start и --end")

    from sam_reader import SamReader
    reader = SamReader(input_file)
    
    if action in (None, "count"):
        click.echo(f"Кол-во выравниваний: {reader.count_records()}")
    if action in (None, "stats"):
        click.echo("Статистика по хромосомам:")
        click.echo(reader.chromosome_stats())
    if action == "region":
        click.echo(f"Выравнивания в регионе {chrom}:{start}-{end}:")
        count = sum(1 for _ in reader.get_region(chrom, start, end))
        click.echo(f"Найдено: {count}")


# КОМАНДА: VCF

@main.command()
@click.option("-i", "--input", "input_file", required=True, help="Входной файл")
@click.option("-a", "--action", type=click.Choice(["count", "stats", "region"]), default=None,
              help="Действие: count (кол-во), stats (статистика), region (поиск в регионе)")
@click.option("--chrom", help="Хромосома для поиска")
@click.option("--start", type=int, help="Начало региона")
@click.option("--end", type=int, help="Конец региона")
def vcf(input_file, action, chrom, start, end):
    """Анализ VCF"""
    if action == "region" and not (chrom and start and end):
        raise click.UsageError("Для действия 'region' необходимо указать --chrom, --start и --end")

    from vcf_reader import VcfReader
    reader = VcfReader(input_file)
    
    if action in (None, "count"):
        click.echo(f"Кол-во вариантов: {reader.count_variants()}")
    if action == "stats":
        click.echo("Статистика по регионам (окно 100000):")
        click.echo(reader.region_count_stats(bin_size=100000))
    if action == "region":
        click.echo(f"Варианты в регионе {chrom}:{start}-{end}:")
        count = sum(1 for _ in reader.get_region(chrom, start, end))
        click.echo(f"Найдено: {count}")


if __name__ == "__main__":
    main()