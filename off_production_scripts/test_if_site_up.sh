#! /usr/bin/sh

#NOTE: Assumes mailutils are already set up (see setup.txt).
#For an SMTP Server, use the Gmail account: darkreactionsproject@gmail.com
EMAIL="darkreactionsproject@lists.haverford.edu"
URL="darkreactions.haverford.edu"

#For more: "http://answers.google.com/answers/threadview/id/276934.html"
function check_email {
  if [ "$1" -ne 0 ] ; then
    #If the error code is not 0 (ie: success)
    echo "Could not access site: $2"
    send_error_email $2 $3
  fi
}


#Check out: "http://stackoverflow.com/questions/8260858/how-to-send-email-from-terminal"
function send_error_email {
  echo "Could not access $1 !" | mail -s "DRP: Error Accessing Site" "$2"
}


curl -s -o "/dev/null" $URL
check_email $? $URL $EMAIL
