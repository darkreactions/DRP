from abc import ABCMeta, abstractmethod

class AbstractRecommender(object):
  __metaclass__ = ABCMeta

  @abstractmethod
  def recommend(self):
    pass


