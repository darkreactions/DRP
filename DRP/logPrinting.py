import sys
import datetime

#Prints an error to stderr for logging purposes.
def print_error(error, details=None):
  sys.stderr.write("{}\n".format(datetime.datetime.now()))
  sys.stderr.write("ERROR: {}\n".format(error))
  if details:
    sys.stderr.write("DETAILS: {}\n".format(details))
  sys.stderr.write("________\n")
  sys.stderr.flush()

def print_log(message):  
  #Note: Using "stdout.write" so that output is flushed programatically.
  sys.stdout.write("{}\n".format(datetime.datetime.now()))
  sys.stdout.write("INFO: {}\n".format(message))
  sys.stdout.write("________\n")
  sys.stdout.flush()
