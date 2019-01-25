GitHub Branches
===============

master
  The current up-to-date working branch, that users download It should
  only contain the latest release and bug fixes.

develop
  It is a fork of the *master* branch with new developments that will be
  merged in the *master* branch for the next release.

gh-pages
  This is an independent branch, containing only the web site of QP2.


# How to make a bug fix

[git-flow](https://nvie.com/posts/a-successful-git-branching-model)
should be used:

[![git-flow]](https://nvie.com/img/git-model@2x.png)

For each bug fix, create a new branch. Then, when the bug is fixed and
the test are all OK, merge (squash) the bug fix commits in both the
*master* and the *develop* branches.


