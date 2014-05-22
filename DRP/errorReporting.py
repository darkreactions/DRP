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
  print "{}\n".format(datetime.datetime.now())
  print "INFO: {}\n".format(message)
  print "________\n"
