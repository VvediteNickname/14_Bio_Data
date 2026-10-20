from abc import ABC, abstractmethod

class Reader(ABC):
    
    def __init__(self, filename):
        self.filename = filename  # Просто сохраняем путь к файлу

    @abstractmethod
    def read(self):
        """Заглушка. Все потомки обязаны написать свою версию этого метода."""
        pass