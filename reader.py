from abc import ABC, abctractmethod

class Reader(ABC):
    
    def __init__(self, filename):
        self.filename = filename  # Просто сохраняем путь к файлу

    @abctractmethod
    def read(self):
        """Заглушка. Все потомки обязаны написать свою версию этого метода."""
        pass