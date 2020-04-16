# Contributing to universc development

Thank you for helping to make this package better. We value all contributions
and rely on your feedback to identify problems and use cases.

## TL;DR

Send your PR! Thanks!

## More Details

You want to contribute? Awesome! Small changes, like fixing typos in
documentation are completely fine and also most welcome. For bigger
changes, we suggest that you open an issue before you start coding, so that
we can maximize the probability that we can successfully merge in your
code.

The goal of this guide is to help you get up and contributing to universc as 
quickly as possible. The guide is divided into two main pieces:
  
  1. Filing a bug report or feature request in an issue.
  
  2. Suggesting a change via a pull request.

Please note that universc is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project,  you agree to abide by its terms.

## Issues

When filing an issue, the most important thing is to include a minimal 
reproducible example so that we can quickly verify the problem, and then figure 
out how to fix it. There are three things you need to include to make your 
example reproducible: required packages, data, code.

1.  **Packages** should be added to the PATH at the top of the script,
    so it's easy to see which ones the example needs. You can see which
    you are using with:

```
echo $PATH
```

```
grep "PATH" $HOME/bash.r  
```

1.  The easiest way to include **data** is to use `head` to show the
    first few rows of the files we are dealing with.

```
head -n 24 *R1_001.fastq
head -n 24 *R2_001.fastq
head -n 10 barcodes.txt
```

You can show compressed files as follows:

```
for file in `ls *R1_001.fastq *R2_001.fastq`
do
    zcat $file | head -n 24
done
```

Then simply copy these into GitHub with code blocks:

\`\`\`

your code here

\`\`\`

1. Please include error messages with the full backtrace of
   commands that failed.

```
tail -n 25  my-id/outs/_log
```

Copy the output from the cellranger logs or the terminal
and paste is as a quote

> \>

> \> your quote here

> \>

1.  Spend a little bit of time ensuring that your **code** is easy for others to
read:
  
  * make sure you've used spaces and your variable names are concise, but
      informative
  
    * use comments to indicate where your problem lies
  
    * do your best to remove everything that is not related to the problem.  
     The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a 
fresh shell session and pasting your script in.

(Unless you've been specifically asked for it, please don't include the
packages that you've installed.)

## Pull requests

To contribute a change to universc, you follow these steps:

1. Create a branch in git and make your changes.
1. Push branch to github and issue pull request (PR).
1. Discuss the pull request.
1. Iterate until either we accept the PR or decide that it's not
a good fit for universc.

Each of these steps are described in more detail below. This might feel 
overwhelming the first time you get set up, but it gets easier with practice. 
If you get stuck at any point, please reach out for help on the [universc-dev](https://groups.google.com/forum/#!forum/universc-dev) mailing list.
                                                                                
If you're not familiar with git or github, please start by reading <http://r-pkgs.had.co.nz/git.html>

<!--
* [ ] Motivate the change in one paragraph, and include it in NEWS.
      In parentheses, reference your github user name and this issue:
      `(@hadley, #1234)`
* [ ] Check pull request only includes relevant changes.
* [ ] Use the [Google shell style](https://google.github.io/styleguide/shellguide.html).
* [ ] Update documentation in `man/launch_universc.sh` and the `help` variable in the `launch_universc.sh` script
* [ ] Add test data and a test job if suggesting a new technology (published data is preferred)
* [ ] Add minimal example for a new feature

--->

Pull requests will be evaluated against a seven point checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivate the
    need for change. Unfortunately I am busy with other projects
    these days, so you need to describe the problem and show
    how your pull request solves it as concisely as possible.

    Also include this motivation in `NEWS` so that when a new release of
    universc comes out it's easy for users to see what's changed. Add your
    item at the top of the file and use markdown for formatting. The
    news item should end with `(@yourGithubUsername, #the_issue_number)`.

1.  __Only related changes__. Before you submit your pull request, please
    check to make sure that you haven't accidentally included any unrelated
    changes. These make it harder to see exactly what's changed, and to
    evaluate any unexpected side effects.

    Each PR corresponds to a git branch, so if you expect to submit
    multiple changes make sure to create multiple branches. If you have
    multiple changes that depend on each other, start with the first one
    and don't submit any others until the first one has been processed.
                                                                              
1.  __Use universc coding style__. Please follow the
[Google shell style](https://google.github.io/styleguide/shellguide.html).
Maintaining a consistent style across the whole code base makes it much easier to
jump into the code. If you're modifying existing universc code that
doesn't follow the style guide, a separate pull request to fix the
style would be greatly appreciated.
                                                                              
1.  If you're adding new parameters or a new function, you'll also need
to document them in the help and manual.
                                                                              
1.  If fixing a bug or adding a new feature please add a
    test to check if it works if test data is available.
                                                                              
1.  If fixing a bug in the output, please give an example.
                                                                              
1.  If you're adding a new feature, please add a short example
    to the appropriate function.

This seems like a lot of work but don't worry if your pull request isn't perfect.
It's a learning process and members of the universc team will be on hand to help you
out. A pull request ("PR") is a process, and unless you've submitted a few in the
past it's unlikely that your pull request will be accepted as is. All PRs require
review and approval from at least one member of the universc development team 
before merge.
                                                                              
Please remember that universc is package used by other people. 
This means that changing any existing functionality could result in
breaking someone's code (or another pipeline). 
Please don't submit pull requests that change existing behaviour. Instead, 
think about how you can add a new feature in a minimally invasive way.

## Making Small Changes

* Please always use the `dev` branch. Choose this branch in your fork. (We
  build the `master` branch from the `dev` branch automatically, to make
  sure that the repo is compatible with the `devtools` R package which uses
  the `master` branch by default.)
* Then look for the file you want to modify.
* Click on the edit symbol (pen) on the upper right corner of the file
  view.
* Make your edits.
* Write a short commit message, less than 65 characters. E.g.  "Fix manual
  page typo" or "Fix degree bug for loops". If needed, elaborate your
  changes below in the "extended description" field.
* Commit your changes.
* Go back to the start page of *your* forked repository. It is at
  `https://github.com/<username>/universc`.
* Click on the green button before the branch name to create a pull
  request.
* Click on "Create pull request".
* Provide a more detailed description if you like. Please also indicate
  that you are fine with licensing your contribution under universc's license
  (see Legal Stuff below).
* Click on "Create pull request".
* That's it! It is probably a good idea to keep your forked repository
  until the change is accepted into universc, in case you need to modify it.
* Now you need to wait for us, unfortunately. Please ping us, if it takes
  long to respond. E.g. a week is considered to be long.
* Once your pull request is accepted, you can delete your forked repository.

## Making More Involved Changes

This is mostly the same as for trivial changes, but you probably want to
edit the sources on your computer, instead of online on Github.

* Open an issue in the issue tracker about the proposed changes.  This is
  not required for smaller things, but I suggest you do it for others. Just
  in case somebody is already working on the same thing, or it is something
  we don't want in universc.
* Fork the repository, and clone it to the machine you'll work on.
* We usually build universc on OSX, so the `dev` branch is usally fine on
  that platform. It might have problems on other systems. If this happens,
  please open an issue and tell us.
* Make sure you work on the `dev` branch.
* Once ready with your changes, build universc, and run the tests. If you use
  the `docker` package, this means running:

  ```sh
  docker build docker build -t universc:dev .
  ```

* Submit your pull request.
* Now you need to wait for us, unfortunately. Please ping us,
  by email or mentioning the maintainer's username on GitHub
  if it takes longer than a week or so to respond. 

## Writing universc Code 

Some tips on writing universc code. In general, look at how things are done,
and try to do them similarly. (Unless you think they are not done well, in
which case please tell us.)

### Code Formatting

Look at the style (indentation, braces, etc.) of some recently committed
bigger change, and try to mimic that. The code style within universc is not
stricly the same, but we want to keep it reasonably similar. If you are 
unsure on this, we can address this when reviewing the Pull Request so
don't worry about it too much.

### Documentation

Please document your new functions by editing the manual.

```
nano man/launch_universc.sh
```

```
nano +$(grep -n "help=" launch_universc.sh | cut -d: -f1) launch_universc.sh
```

### Test Cases

Unless you change something trivial, please consider adding test cases.
This is important! See the files in the `test` directory for
examples.

If you are unsure how to do this please submit the PR anyway and
we can discuss it. If you are able to describe the use case
where this feature would be useful, please do so as best you can.

### Ask Us!

In general, if you are not sure about something, please ask! You can
open an issue on Github or contact the maintainer <tom.kelly@riken.jp>.
We to answer publicly so that others can learn from it, too. There
are silly questions, if you're having trouble others probably are too.

## Legal Stuff

This is a pain to deal with, but we can't avoid it, unfortunately.  So,
universc is licensed under the "General Public License (GPL) version 3, or
later". If your contribution is bigger than a typo fix, then please
indicate that you are fine with releasing your code/text under these
licenses.  E.g. adding a sentence that reads as "I'm fine with GPL 3"
is perfectly enough.
