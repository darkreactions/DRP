from AbstractSplitter import AbstractSplitter
import random

class Splitter(AbstractSplitter):
  def __init__(self, namingStub):
    super(Splitter, self).__init__(namingStub)

  def split(self, reactions):
    splits = [ (self.package(reactions), self.package([])) ]

    return splits
