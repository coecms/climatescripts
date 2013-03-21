climatescripts
==============

Useful functions from researchers within the CoE and CCRC

Contributing
------------

1. First go to https://github.com/ and register
2. Then go the COE repository: https://github.com/coecms/climatescripts
3. On the top right hand side, click on the 'Fork' button, and you should see something like:
https://github.com/"Your_Github_username"/climatescripts.git
4. On your local machine, create a directory where you will keep a copy of the repository
5. Follow the instructions from github on setting ssh-key transfer:
https://help.github.com/articles/generating-ssh-keys
6. Once done, cd to that directory, and Checkout the repository with git clone:
    
       git clone git@github.com:"Your_Github_username"/climatescripts.git

7. Edit existing files or add your own
8. Run git status - This tells you what files are being tracked, have changed etc.
9. If you want to add new script file, use git add "file_name".
9.  Commit: Commit your changes and push them to your github account
 
        git commit -a -m "A useful message"
    
        git push

6. Pull request: Go to  https://github.com/coecms/climatescripts, Press the pull request button on github (you will need to add a mesage). This sends a message to the admin person that you have added stuff, and that person decides whether to merge the code into the repository

Updating your local repository
------------------------------

To get the latest changes from the master repository you'll first need to tell git its address

     git remote add coecms git://github.com/coecms/climatescripts.git

then update by running 

    git pull coecms master

Some Considerations
-------------------

Please comment/document and indent your code - In Python, you do not get a choice, but not all other languages
Please review other people's code whenever possible
Very good basic GIT tutorials can be found at: http://git-scm.com/book/en/Getting-Started-Git-Basics
In a netshell: git status; git add; git commit; git push; pull request
It is usually better to commit incrementally, rather than once-off with lots of changes to the code 
