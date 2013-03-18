climatescripts
==============

Useful functions from researchers within the CoE and CCRC

Contributing
------------

1. Fork: Click the 'Fork' button on github.com/coecms/climatescripts
2. Checkout: Checkout the repository with git clone, e.g.

    git clone git@github.com:ScottWales/climatescripts

3. Branch: Create a branch for new functionality

    git checkout -b 'fieldsfile'

4. Edit: Add in your changes
5. Commit: Commit your changes and push them to your github account

    git commit -a -m "Added new fieldsfile functions"
    
    git push origin fieldsfile

6. Pull request: Press the pull request button on github and select the branch to send

Updating your local repository
------------------------------

To get the latest changes from the master repository you'll first need to tell git its address

     git remote add coecms git://github.com/coecms/climatescripts.git

then update by running 

    git pull coecms master
