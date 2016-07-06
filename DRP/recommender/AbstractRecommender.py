from abc import ABCMeta, abstractmethod

class AbstractRecommender(object, metaclass=ABCMeta):
  @abstractmethod
  def recommend(self):
    pass


