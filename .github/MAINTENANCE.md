# Maintenance

The study package tarball is used by the OmicNavigator [quick start
instructions][quick] to demo the app's features. Whenever you make a change to
the study package, make sure to regenerate the tarball using the code below.
It's important to remove the version number for the tarball file name so that
the quick start instructions don't have to be updated every time the RNAseq123
study package is updated.

[quick]: https://github.com/abbvie-external/OmicNavigator#quick-start

```R
library("OmicNavigator")
x <- importStudy("RNAseq123")
tarball <- exportStudy(x, type = "tarball")
file.rename(tarball, "ONstudyRNAseq123.tar.gz")
```

Then create a new tag/release and upload `ONstudyRNAseq123.tar.gz` as a release
asset. The latest tarball will always be available to download as
https://github.com/abbvie-external/OmicNavigatorExample/releases/latest/download/ONstudyRNAseq123.tar.gz
([GitHub docs][latest]).

[latest]: https://docs.github.com/en/github/administering-a-repository/releasing-projects-on-github/linking-to-releases
