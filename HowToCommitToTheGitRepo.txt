#  These are instructions to set up BitBucket to use an SSH key
#  so that we can have multiple accounts pushing to the same repo.
#  Below, replace "cfalk" with whatever your username is.

# # # # Server Side # # # # #
#Start up the ssh-agent in bash.
exec ssh-agent bash

#Add the RSA key you got to a file. You want to have a passphrase.
ssh-keygen

#Yield the contents of the file it created. Go ahead and copy them.
cat /home/cfalk/.ssh/id_rsa.pub


# # # # BitBucket Side # # # #
#Now go to your personal account settings on BitBucket.

#Click the "SSH keys" setting to the left.

#Paste the id_rsa.pub contents in and click "Add Key".

#You should be able to Git as you please!


#For more thoroughness, see: https://confluence.atlassian.com/pages/viewpage.action;jsessionid=D25198F62DD02A07968C6B4AC6F35AEB?pageId=270827678 .

