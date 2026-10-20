import os
import matplotlib.pyplot as plt
from reader import Reader 

class SequenceReader(Reader):
    
    def validate_sequence(self, sequence):  #проверка ошибок: смотрим, чтобы в ДНК были только буквы A, T, G, C
        sequence = sequence.upper() # перевод строки в большие буквы
        for letter in sequence: # проверяем каждую букву
            if letter not in ['A', 'T', 'G', 'C', 'N']:
                return False  # если нашли не ту букву — ошибка
        return True

    def get_sequences_count(self):  #Считаем, сколько всего штук ДНК в файле
        if not os.path.exists(self.filename):           #файла нет — ошибка
            print("Ошибка: Файл не найден!")
            return 0
            
        count = 0
        # перебираем всё, что выдало
        for item in self.read():
            count = count + 1  # нашли 1 ДНК — прибавили 1
        return count

    def get_average_length(self):   #считаем среднюю длину ДНК в файле
        if not os.path.exists(self.filename):
            return 0
            
        total_length = 0
        count = 0
        
        for item in self.read():
            dna_string = item[1] # строка с ДНК идет второй в результате
            total_length = total_length + len(dna_string) # считаем длину этой строки и прибавляем к общей сумме
            count = count + 1
            
        # если файл пустой, чтобы не делить на ноль
        if count == 0:
            return 0
            
        return total_length / count # делим общую длину на количество штук

    def plot_quality(self): #график качества (только FASTQ) для ДНК любой длины
        if not os.path.exists(self.filename):
            print("Ошибка: Файл не найден!")
            return

        all_positions = {} #создаем словарь

        for item in self.read():
            # если в результате <3 элементов, значит это FASTA (там нет качества), выходим (имя, днк, качество)
            if len(item) < 3:
                print("График качества можно строить только для FASTQ файлов!")
                return
                
            quality_string = item[2]  # строка со значками качества
            
            for index in range(len(quality_string)): # идем по каждому значку в строке качества
                char = quality_string[index]
                number_score = ord(char) - 33 # переводим значок в число по биологической формуле. 
                #ord переводит символ в порядк номер в памяти компьютера !=33
                
                # если такой позиции еще не было в словаре, создаем ее
                if index not in all_positions:
                    all_positions[index] = []
                
                # добавляем оценку качества для этой позиции
                all_positions[index].append(number_score)

        #Проверяем, есть ли данные
        if not all_positions:
            print("Файл пуст или не содержит данных о качестве.")
            return

        # считаем среднее число для каждой позиции
        x_points = []  # номера позиций (0, 1, 2...)
        y_points = []  # средние оценки качества
        
        # сортируем позиции по порядку (0, 1, 2, 3...)
        sorted_positions = sorted(all_positions.keys())
        
        for index in sorted_positions:
            list_of_scores = all_positions[index]
            # считаем среднее для этой позиции, сумму делим на кол-во чисел
            average = sum(list_of_scores) / len(list_of_scores)
            x_points.append(index)
            y_points.append(average)

        # рисуем график
        plt.figure()
        plt.plot(x_points, y_points, color="blue")
        plt.title("Качество FASTQ")
        plt.xlabel("Позиция")
        plt.ylabel("Качество")
        
        # сохраняем график
        plt.savefig("quality_plot.png")
        plt.close()
        print("График успешно сохранен в файл quality_plot.png")

class FastaReader(SequenceReader):
    
    def read(self): #читаем файлы FASTA
        if not os.path.exists(self.filename):
            return

        name = None
        dna_letters = []

        with open(self.filename, "r") as f:
            for line in f:
                line = line.strip()  # убираем пробелы в конце строки
                
                if line == "":  # если строка пустая, пропускаем её
                    continue

                if line.startswith(">"):  # если строка начинается с '>', это имя новой ДНК
                    if name is not None:
                        full_dna = "".join(dna_letters)  # склеиваем все буквы в один текст
                        # проверка ошибок:
                        if self.validate_sequence(full_dna) == False:
                            print("Внимание! В ДНК нашли другие буквы!")
                        yield name, full_dna  # отдаем имя и ДНК
                    
                    name = line[1:]  # запоминаем имя (без'>')
                    dna_letters = []  # очищаем список букв для следующей ДНК
                else:
                    # если строка обычная, это буквы ДНК, добавляем их в список
                    dna_letters.append(line)

            # отдаем последнюю ДНК из файла
            if name is not None:
                full_dna = "".join(dna_letters)
                if self.validate_sequence(full_dna) == False:
                    print("Внимание! В ДНК нашли другие буквы!")
                yield name, full_dna


class FastqReader(SequenceReader):
    
    def read(self):  #читаем файлы FASTQ блоками по 4 строки
        if not os.path.exists(self.filename):
            return

        with open(self.filename, "r") as f:
            while True:
                line1 = f.readline()
                line2 = f.readline()
                line3 = f.readline()
                line4 = f.readline()

                # если первая строка пустая, значит файл закончился, выходим из цикла
                if line1 == "":
                    break
                # проверка ошибок: если файл оборвался на середине блока
                if line2 == "" or line3 == "" or line4 == "":
                    print("Ошибка: Файл сломан, в нем не хватает строк!")
                    break

                name = line1.strip()[1:]      # имя (убираем '@')
                dna = line2.strip()           # буквы ДНК
                quality = line4.strip()       # значки качества

                # Проверка ошибок:
                if self.validate_sequence(dna) == False:
                    print("Внимание! В FASTQ нашли другие буквы!")
                if len(dna) != len(quality):
                    print("Ошибка: Кол-во букв и значков качества не совпадает!")
                yield name, dna, quality
